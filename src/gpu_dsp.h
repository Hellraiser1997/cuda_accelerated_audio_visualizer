#ifndef GPU_DSP_H
#define GPU_DSP_H

#include <stdint.h>
#include <stdbool.h>

// Initialize any GPU resources needed
bool init_gpu_dsp();

// Process a chunk of audio samples on the GPU and get back logarithmic frequency magnitudes
// audio_data: 16-bit PCM samples
// num_frames: number of frames to process (e.g., 2048)
// num_channels: 1 for Mono, 2 for Stereo
// out_magnitudes: array to store the resulting frequency magnitudes (interleaved L, R if stereo)
// out_num_bins: number of frequency bins returned per channel
bool process_audio_chunk(const int16_t* audio_data, uint32_t num_frames, uint16_t num_channels, float* out_magnitudes, uint32_t* out_num_bins);

// Process a chunk of audio samples on the GPU and get back the 12 musical pitch classes (Chromagram)
// out_chroma: array to store the 12 pitch magnitudes (interleaved L, R if stereo)
// out_num_pitches: always returns 12
bool compute_chromagram(const int16_t* audio_data, uint32_t num_frames, uint16_t num_channels, float* out_chroma, uint32_t* out_num_pitches);

// Cleanup GPU resources
void cleanup_gpu_dsp();

#endif // GPU_DSP_H
