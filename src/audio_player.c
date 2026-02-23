#include "audio_player.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

#ifdef USE_PULSEAUDIO
#include <pulse/simple.h>
#include <pulse/error.h>
pa_simple *pa_handle = NULL;
#else
#include <alsa/asoundlib.h>
snd_pcm_t *pcm_handle = NULL;
// snd_pcm_t *pcm_handle = NULL; // This is now part of AudioPlayer
#endif

pthread_t playback_thread;

void* playback_loop(void* arg) {
    AudioPlayer* player = (AudioPlayer*)arg;
    AudioData* data = player->audio_data;
    uint32_t total_frames = data->num_samples; // which is actually frames based on media_reader output
    
    // Frames to process each iteration
    uint32_t frames_to_write_chunk = 1024;
    
    while (player->is_playing && player->current_sample_index < total_frames) {
        uint32_t remaining = total_frames - player->current_sample_index;
        uint32_t frames_to_write = frames_to_write_chunk;
        if (remaining < frames_to_write) frames_to_write = remaining;
        
        size_t bytes_to_write = frames_to_write * data->header.num_channels * 2; // 2 bytes per sample (S16)
        size_t array_offset = player->current_sample_index * data->header.num_channels;
        
#ifdef USE_PULSEAUDIO
        pa_simple* pa_handle = (pa_simple*)player->backend_handle;
        if (!pa_handle) {
            uint32_t sleep_us = (uint32_t)(((float)frames_to_write / TARGET_SAMPLE_RATE) * 1000000.0f);
            usleep(sleep_us);
            player->current_sample_index += frames_to_write;
            continue;
        }

        int pa_error;
        if (pa_simple_write(pa_handle, &data->audio_data[array_offset], bytes_to_write, &pa_error) < 0) {
            // PulseAudio failed, or stream closed midway
            fprintf(stderr, "pa_simple_write() failed: %s\n", pa_strerror(pa_error));
            break;
        } else {
            player->current_sample_index += frames_to_write;
        }
#else
        snd_pcm_t* pcm_handle = (snd_pcm_t*)player->backend_handle;
        if (!pcm_handle) {
            uint32_t sleep_us = (uint32_t)(((float)frames_to_write / TARGET_SAMPLE_RATE) * 1000000.0f);
            usleep(sleep_us);
            player->current_sample_index += frames_to_write;
            continue;
        }

        // Write interleaved data to ALSA
        int err = snd_pcm_writei(pcm_handle, &data->audio_data[array_offset], frames_to_write);
        
        if (err == -EPIPE) {
            // Buffer underrun
            fprintf(stderr, "ALSA underrun occurred\n");
            snd_pcm_prepare(pcm_handle);
        } else if (err < 0) {
            fprintf(stderr, "Error writing to PCM device: %s\n", snd_strerror(err));
            player->is_playing = false;
        } else {
            player->current_sample_index += err;
        }
#endif
    }
    
    player->is_playing = false;

#ifdef USE_PULSEAUDIO
    // pa_simple_drain is handled implicitly or better handled by the main thread cleanup
    // doing it here while the main thread might be calling pa_simple_free causes an abort
#else
    if (player->backend_handle) {
        snd_pcm_drain((snd_pcm_t*)player->backend_handle);
    }
#endif

    return NULL;
}

bool init_audio_player(AudioPlayer* player, AudioData* data) {
    if (!player || !data) return false;
    
    player->is_playing = false;
    player->backend_handle = NULL;
    player->audio_data = data;
    player->current_sample_index = 0;

#ifdef USE_PULSEAUDIO
    // PulseAudio setup
    static const pa_sample_spec ss = {
        .format = PA_SAMPLE_S16LE,
        .rate = TARGET_SAMPLE_RATE, // standardized format
        .channels = TARGET_CHANNELS
    };
    int pa_error;
    // Open a simple PulseAudio stream for playback
    pa_simple* pa_handle_local = pa_simple_new(NULL,               // Use the default server.
                              "CudaMusicPlayer",  // Our application's name.
                              PA_STREAM_PLAYBACK,
                              NULL,               // Use the default device.
                              "Music",            // Description of our stream.
                              &ss,                // Our sample format.
                              NULL,               // Use default channel map
                              NULL,               // Use default buffering attributes.
                              &pa_error);         // Ignore error code.

    if (!pa_handle_local) {
        fprintf(stderr, "Warning: pa_simple_new() failed: %s. Falling back to visual-only mode.\n", pa_strerror(pa_error));
        // We still return true to allow the visualizer to run in mocked time
        return true; 
    }
    player->backend_handle = (void*)pa_handle_local;
    printf("PulseAudio backend initialized successfully.\n");
#else
    int err;
    snd_pcm_t* pcm_handle_local = NULL;
    // Open ALSA device "default"
    if ((err = snd_pcm_open(&pcm_handle_local, "default", SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
        fprintf(stderr, "Warning: Cannot open ALSA device default (%s). Falling back to visual-only mode.\n", snd_strerror(err));
        player->backend_handle = NULL;
        return true; 
    }
    player->backend_handle = (void*)pcm_handle_local;

    // Allocate parameters object and fill it with default values
    snd_pcm_hw_params_t *params;
    snd_pcm_hw_params_alloca(&params);
    snd_pcm_hw_params_any(pcm_handle_local, params);

    // Interleaved mode
    snd_pcm_hw_params_set_access(pcm_handle_local, params, SND_PCM_ACCESS_RW_INTERLEAVED);

    // Signed 16-bit little-endian format
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_hw_params_set_format(pcm_handle_local, params, format);

    // Channels
    snd_pcm_hw_params_set_channels(pcm_handle_local, params, data->header.num_channels);

    // Sample rate
    unsigned int rate = data->header.sample_rate;
    snd_pcm_hw_params_set_rate_near(pcm_handle_local, params, &rate, 0);

    // Write parameters
    if ((err = snd_pcm_hw_params(pcm_handle_local, params)) < 0) {
        fprintf(stderr, "Cannot set ALSA parameters (%s)\n", snd_strerror(err));
        return false;
    }
    printf("ALSA backend initialized successfully.\n");
#endif

    return true;
}

void start_playback(AudioPlayer* player) {
    if (!player || player->is_playing) return;
    
    player->is_playing = true;
    
    // Spawn background thread to feed Audio backend
    if (pthread_create(&playback_thread, NULL, playback_loop, player) != 0) {
        fprintf(stderr, "Failed to create playback thread\n");
        player->is_playing = false;
    }
}

void cleanup_audio_player(AudioPlayer* player) {
    if (player->is_playing) {
        player->is_playing = false;
        pthread_join(playback_thread, NULL);
    }
    
#ifdef USE_PULSEAUDIO
    if (pa_handle) {
        pa_simple_free(pa_handle);
        pa_handle = NULL;
    }
#else
    if (pcm_handle) {
        snd_pcm_close(pcm_handle);
        pcm_handle = NULL;
    }
#endif
}
