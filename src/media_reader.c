#include "media_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libswresample/swresample.h>
#include <libavutil/opt.h>

bool load_media(const char* filepath, AudioData* out_data) {
    if (!filepath || !out_data) return false;

    AVFormatContext* format_context = NULL;
    AVCodecContext* codec_context = NULL;
    const AVCodec* codec = NULL;
    SwrContext* swr_context = NULL;
    AVFrame* frame = NULL;
    AVPacket* packet = NULL;
    
    int audio_stream_index = -1;
    bool success = false;

    // 1. Open the file and read header information
    if (avformat_open_input(&format_context, filepath, NULL, NULL) != 0) {
        fprintf(stderr, "FFmpeg: Could not open file %s\n", filepath);
        goto cleanup;
    }

    if (avformat_find_stream_info(format_context, NULL) < 0) {
        fprintf(stderr, "FFmpeg: Could not find stream info\n");
        goto cleanup;
    }

    // 2. Find the best audio stream
    audio_stream_index = av_find_best_stream(format_context, AVMEDIA_TYPE_AUDIO, -1, -1, &codec, 0);
    if (audio_stream_index < 0) {
        fprintf(stderr, "FFmpeg: Could not find any audio stream in the file\n");
        goto cleanup;
    }

    // 3. Allocate and initialize the codec context
    codec_context = avcodec_alloc_context3(codec);
    if (!codec_context) {
        fprintf(stderr, "FFmpeg: Failed to allocate codec context\n");
        goto cleanup;
    }

    if (avcodec_parameters_to_context(codec_context, format_context->streams[audio_stream_index]->codecpar) < 0) {
        fprintf(stderr, "FFmpeg: Failed to copy codec parameters to context\n");
        goto cleanup;
    }

    if (avcodec_open2(codec_context, codec, NULL) < 0) {
        fprintf(stderr, "FFmpeg: Failed to open codec\n");
        goto cleanup;
    }

    // Ensure channel layout is valid on the input
    if (codec_context->ch_layout.nb_channels == 0) {
       fprintf(stderr, "FFmpeg: Unknown input channel layout\n");
       goto cleanup;
    }

    // 4. Set up the SwrContext for forced resampling into 44.1kHz Stereo 16-bit Int
    AVChannelLayout out_ch_layout;
    av_channel_layout_default(&out_ch_layout, TARGET_CHANNELS);
    
    enum AVSampleFormat out_sample_fmt = AV_SAMPLE_FMT_S16; // 16-bit signed interleaved
    int out_sample_rate = TARGET_SAMPLE_RATE;

    if (swr_alloc_set_opts2(&swr_context,
                            &out_ch_layout, out_sample_fmt, out_sample_rate,
                            &codec_context->ch_layout, codec_context->sample_fmt, codec_context->sample_rate,
                            0, NULL) < 0) {
        fprintf(stderr, "FFmpeg: Failed to allocate SWResample context\n");
        goto cleanup;
    }

    if (swr_init(swr_context) < 0) {
        fprintf(stderr, "FFmpeg: Failed to initialize SWResample context\n");
        goto cleanup;
    }

    // 5. Allocate buffers
    frame = av_frame_alloc();
    packet = av_packet_alloc();
    if (!frame || !packet) {
        fprintf(stderr, "FFmpeg: Failed to allocate frame or packet\n");
        goto cleanup;
    }

    // Estimate output memory (Duration * Sample Rate * Channels * 2 bytes)
    // We add 10% padding because VBR files can sometimes miscalculate exact duration
    int64_t est_duration = format_context->duration != AV_NOPTS_VALUE ? format_context->duration : 0;
    int64_t est_total_samples = (est_duration * out_sample_rate) / AV_TIME_BASE;
    if (est_total_samples == 0) est_total_samples = out_sample_rate * 300; // default 5 mins if unknown
    
    size_t allocated_bytes = (size_t)(est_total_samples * 1.1) * TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE;
    int16_t* pcm_buffer = (int16_t*)malloc(allocated_bytes);
    if (!pcm_buffer) {
        fprintf(stderr, "FFmpeg: Failed to allocate PCM output buffer\n");
        goto cleanup;
    }

    size_t current_pcm_offset_bytes = 0;
    
    // 6. Decode Loop
    while (av_read_frame(format_context, packet) >= 0) {
        if (packet->stream_index == audio_stream_index) {
            
            // Send packet to decoder
            int ret = avcodec_send_packet(codec_context, packet);
            if (ret < 0) {
                av_packet_unref(packet);
                break;
            }

            // Receive all decoded frames from this packet
            while (ret >= 0) {
                ret = avcodec_receive_frame(codec_context, frame);
                if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF) {
                    break;
                } else if (ret < 0) {
                    fprintf(stderr, "FFmpeg: Error during decoding\n");
                    goto decode_error;
                }

                // Resample the frame
                // Find out how many output samples will be generated
                int out_samples_per_channel = av_rescale_rnd(
                    swr_get_delay(swr_context, codec_context->sample_rate) + frame->nb_samples,
                    out_sample_rate, codec_context->sample_rate, AV_ROUND_UP);
                
                int out_buffer_size_needed = out_samples_per_channel * TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE;
                
                // Reallocate if our estimated buffer was too small
                if (current_pcm_offset_bytes + out_buffer_size_needed > allocated_bytes) {
                    allocated_bytes = allocated_bytes * 2 + out_buffer_size_needed; // Double it + ensure enough space
                    int16_t* new_buffer = (int16_t*)realloc(pcm_buffer, allocated_bytes);
                    if (!new_buffer) {
                        fprintf(stderr, "FFmpeg: Reallocation failed\n");
                        goto decode_error;
                    }
                    pcm_buffer = new_buffer;
                }

                // Actually run the conversion
                uint8_t* out_ptr = (uint8_t*)pcm_buffer + current_pcm_offset_bytes;
                int converted_samples_per_channel = swr_convert(swr_context,
                                                    &out_ptr, out_samples_per_channel,
                                                    (const uint8_t**)frame->data, frame->nb_samples);
                                                    
                if (converted_samples_per_channel < 0) {
                    fprintf(stderr, "FFmpeg: Error while resampling\n");
                    goto decode_error;
                }

                current_pcm_offset_bytes += converted_samples_per_channel * TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE;
                av_frame_unref(frame);
            }
        }
        av_packet_unref(packet);
    }
    
    // Process any remaining delayed samples in the resampler
    int flush_samples = swr_get_out_samples(swr_context, 0);
    if (flush_samples > 0) {
        int out_buffer_size_needed = flush_samples * TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE;
        if (current_pcm_offset_bytes + out_buffer_size_needed > allocated_bytes) {
            allocated_bytes += out_buffer_size_needed;
            pcm_buffer = (int16_t*)realloc(pcm_buffer, allocated_bytes);
        }
        uint8_t* out_ptr = (uint8_t*)pcm_buffer + current_pcm_offset_bytes;
        int converted_samples = swr_convert(swr_context, &out_ptr, flush_samples, NULL, 0);
        if (converted_samples > 0) {
             current_pcm_offset_bytes += converted_samples * TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE;
        }
    }

    // 7. Finalize output
    out_data->header.sample_rate = TARGET_SAMPLE_RATE;
    out_data->header.num_channels = TARGET_CHANNELS;
    out_data->header.bits_per_sample = TARGET_FORMAT_BYTES_PER_SAMPLE * 8;
    out_data->audio_data = pcm_buffer;
    out_data->num_samples = current_pcm_offset_bytes / (TARGET_CHANNELS * TARGET_FORMAT_BYTES_PER_SAMPLE);

    success = true;
    goto cleanup;

decode_error:
    free(pcm_buffer);
    success = false;

cleanup:
    if (packet) av_packet_free(&packet);
    if (frame) av_frame_free(&frame);
    if (swr_context) swr_free(&swr_context);
    if (codec_context) avcodec_free_context(&codec_context);
    if (format_context) avformat_close_input(&format_context);
    
    return success;
}

void free_media(AudioData* data) {
    if (data && data->audio_data) {
        free(data->audio_data);
        data->audio_data = NULL;
        data->num_samples = 0;
    }
}
