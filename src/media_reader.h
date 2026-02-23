#ifndef MEDIA_READER_H
#define MEDIA_READER_H

#include <stdint.h>
#include <stdbool.h>

// Standardized Output Format:
// Sample Rate: 44100 Hz
// Channels: 2 (Stereo)
// Format: 16-bit signed integer (Interleaved L R L R)
#define TARGET_SAMPLE_RATE 44100
#define TARGET_CHANNELS 2
#define TARGET_FORMAT_BYTES_PER_SAMPLE 2 // 16-bit

typedef struct {
    uint32_t sample_rate;
    uint16_t num_channels;
    uint16_t bits_per_sample;
} AudioHeader;

typedef struct {
    AudioHeader header;
    int16_t* audio_data; // Pre-allocated and filled Resampled 16-bit PCM buffer
    uint32_t num_samples; // Number of *frames* (1 frame = 1 L sample + 1 R sample)
} AudioData;

// Load, decode, and resample any audio file into standardized 44.1kHz Stereo 16-bit PCM
bool load_media(const char* filepath, AudioData* out_data);

// Free the allocated audio buffer
void free_media(AudioData* data);

#endif // MEDIA_READER_H
