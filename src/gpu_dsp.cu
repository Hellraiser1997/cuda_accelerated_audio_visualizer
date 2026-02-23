#include "gpu_dsp.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

// Logarithmic Frequency Band mapping kernel
// Instead of returning thousands of linear bins, this groups them logarithmically
// which perfectly aligns with human hearing (more detail in bass, groupings in treble).
__global__ void compute_stereo_dft_log_bands(const int16_t* d_audio, float* d_magnitudes, int num_frames, int num_channels, int display_bands, float sample_rate) {
    extern __shared__ int16_t s_audio[];
    
    // Cooperative memory load
    int total_samples = num_frames * num_channels;
    int tid = threadIdx.x;
    int stride = blockDim.x;
    
    for (int i = tid; i < total_samples; i += stride) {
        s_audio[i] = d_audio[i];
    }
    
    __syncthreads(); // Wait until the entire audio chunk is built in shared cache

    int band = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (band < display_bands) {
        
        // Logarithmic frequency bounds for this specific "Band" (e.g., 20Hz to 16,000Hz)
        float min_freq = 20.0f;
        float max_freq = 16000.0f; // Constrain to active musical spectrum rather than full Nyquist
        
        // Calculate the frequency range this band represents using log10
        float log_min = log10f(min_freq);
        float log_max = log10f(max_freq);
        float log_step = (log_max - log_min) / display_bands;
        
        float freq_start = powf(10.0f, log_min + band * log_step);
        float freq_end = powf(10.0f, log_min + (band + 1) * log_step);
        
        // Convert frequencies to discrete FFT bin indices
        int bin_start = (int)((freq_start / sample_rate) * num_frames);
        int bin_end = (int)((freq_end / sample_rate) * num_frames);
        
        if (bin_end <= bin_start) bin_end = bin_start + 1; // Ensure at least 1 bin
        
        // We calculate magnitudes for Left (0) and Right (1) channels
        for (int ch = 0; ch < num_channels; ++ch) {
            float band_magnitude_sum = 0.0f;
            int bins_accumulated = 0;
            
            for (int k = bin_start; k < bin_end; ++k) {
                float real_sum = 0.0f;
                float imag_sum = 0.0f;
                
                // Calculate DFT for bin 'k'
                for (int n = 0; n < num_frames; ++n) {
                    float angle = 2.0f * 3.1415926535f * k * n / num_frames;
                    
                    // Hardware intrinsic math for simultaneous sin/cos execution
                    float s_val, c_val;
                    __sincosf(angle, &s_val, &c_val);
                    
                    // Extract the sample from ultra-fast shared memory cache
                    int sample_idx = n * num_channels + ch;
                    float sample = s_audio[sample_idx] / 32768.0f; 
                    
                    // Apply Hann window
                    float window = 0.5f * (1.0f - cosf(2.0f * 3.1415926535f * n / (num_frames - 1)));
                    sample *= window;
                    
                    real_sum += sample * c_val;
                    imag_sum -= sample * s_val;
                }
                
                float magnitude = sqrtf(real_sum * real_sum + imag_sum * imag_sum);
                band_magnitude_sum += magnitude;
                bins_accumulated++;
            }
            
            // Average the accumulated bins for this band
            float avg_magnitude = (bins_accumulated > 0) ? (band_magnitude_sum / bins_accumulated) : 0;
            
            // Write to output array (Interleave L, R magnitudes for the CPU)
            // Format: [Band0_L, Band0_R, Band1_L, Band1_R, ...]
            d_magnitudes[band * num_channels + ch] = avg_magnitude;
        }
    }
}

// Global pointers for GPU memory
int16_t* d_audio_buffer_bands = NULL;
float* d_magnitude_buffer_bands = NULL;
uint32_t current_buffer_frames_bands = 0;
uint16_t current_channels_bands = 0;

int16_t* d_audio_buffer_chroma = NULL;
float* d_chroma_buffer = NULL;
uint32_t current_buffer_frames_chroma = 0;
uint16_t current_channels_chroma = 0;

bool init_gpu_dsp() {
    int deviceCount;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    if (error_id != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error_id));
        return false;
    }
    return true;
}

// display_bands controls how many physical bars the visualizer will have.
// We can hardcode this to 80 for now since the GUI expects 80.
#define DISPLAY_BANDS 512

bool process_audio_chunk(const int16_t* audio_data, uint32_t num_frames, uint16_t num_channels, float* out_magnitudes, uint32_t* out_num_bands) {
    if (!audio_data || !out_magnitudes || num_frames == 0 || num_channels < 1 || num_channels > 2) return false;

    if (out_num_bands) {
        *out_num_bands = DISPLAY_BANDS;
    }

    // Allocate or reallocate Device memory if the chunk size/channels changes
    if (current_buffer_frames_bands != num_frames || current_channels_bands != num_channels) {
        if (d_audio_buffer_bands) cudaFree(d_audio_buffer_bands);
        if (d_magnitude_buffer_bands) cudaFree(d_magnitude_buffer_bands);
        
        cudaMalloc((void**)&d_audio_buffer_bands, num_frames * num_channels * sizeof(int16_t));
        cudaMalloc((void**)&d_magnitude_buffer_bands, DISPLAY_BANDS * num_channels * sizeof(float));
        
        current_buffer_frames_bands = num_frames;
        current_channels_bands = num_channels;
    }

    // 1. Copy interleaved audio chunk from Host to Device
    cudaMemcpy(d_audio_buffer_bands, audio_data, num_frames * num_channels * sizeof(int16_t), cudaMemcpyHostToDevice);

    // 2. Launch Kernel (1 Thread per output Band)
    int blockSize = 256;
    int numBlocks = (DISPLAY_BANDS + blockSize - 1) / blockSize;
    
    // Sample rate is fixed to 44100 for our math currently based on WAV parser expectations
    size_t shared_mem_size = num_frames * num_channels * sizeof(int16_t);
    compute_stereo_dft_log_bands<<<numBlocks, blockSize, shared_mem_size>>>(d_audio_buffer_bands, d_magnitude_buffer_bands, num_frames, num_channels, DISPLAY_BANDS, 44100.0f);
    
    // Wait for GPU to finish
    cudaDeviceSynchronize();

    // 3. Copy Results (magnitudes) from Device to Host
    // Will be [Band0_L, Band0_R, Band1_L, ...]
    cudaMemcpy(out_magnitudes, d_magnitude_buffer_bands, DISPLAY_BANDS * num_channels * sizeof(float), cudaMemcpyDeviceToHost);

    return true;
}

// -----------------------------------------------------------------------------
// CHROMAGRAM (Musical Pitch) KERNEL LOGIC
// -----------------------------------------------------------------------------

// Computes the 12 chromatic pitches across multiple octaves
__global__ void compute_stereo_chromagram(const int16_t* d_audio, float* d_chroma, int num_frames, int num_channels, float sample_rate) {
    extern __shared__ int16_t s_audio_chroma[];
    
    // Cooperative memory load
    int total_samples = num_frames * num_channels;
    int tid = threadIdx.x;
    int stride = blockDim.x;
    
    for (int i = tid; i < total_samples; i += stride) {
        s_audio_chroma[i] = d_audio[i];
    }
    
    __syncthreads();

    // 12 pitch classes (0 = C, 1 = C#, 2 = D ... 11 = B)
    int pitch_class = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (pitch_class < 12) {
        
        // We will scan across 5 octaves starting from C2 (~65.41 Hz) up to C7 (~2093 Hz)
        // A4 = 440Hz -> C2  is 65.40639
        // n = 12 * log2(f / 440) -> f = 440 * 2^(n/12)
        
        for (int ch = 0; ch < num_channels; ++ch) {
            float total_pitch_magnitude = 0.0f;
            
            // Octave 0 starts at C2 (n = -33 relative to A4)
            for (int octave = 0; octave < 6; ++octave) {
                // Calculate the exact target frequency for this specific pitch_class in this octave
                int n_relative_to_A4 = (pitch_class - 9) + (octave - 3) * 12; // C is -9 from A
                float target_freq = 440.0f * powf(2.0f, (float)n_relative_to_A4 / 12.0f);
                
                // Which FFT bin does this perfectly map to?
                int target_bin = (int)((target_freq / sample_rate) * num_frames);
                
                // Allow a small window of bins around the exact target to capture slight out-of-tune harmonics
                int bin_start = target_bin - 2;
                int bin_end = target_bin + 2;
                if (bin_start < 1) bin_start = 1;
                
                float max_bin_magnitude = 0.0f;
                
                for (int k = bin_start; k <= bin_end; ++k) {
                    float real_sum = 0.0f;
                    float imag_sum = 0.0f;
                    
                    // Direct DFT for this specific frequency bin
                    for (int n = 0; n < num_frames; ++n) {
                        float angle = 2.0f * 3.1415926535f * k * n / num_frames;
                        
                        float s_val, c_val;
                        __sincosf(angle, &s_val, &c_val);
                        
                        float sample = s_audio_chroma[n * num_channels + ch] / 32768.0f; 
                        
                        // Hann window
                        float window = 0.5f * (1.0f - cosf(2.0f * 3.1415926535f * n / (num_frames - 1)));
                        sample *= window;
                        
                        real_sum += sample * c_val;
                        imag_sum -= sample * s_val;
                    }
                    
                    float magnitude = sqrtf(real_sum * real_sum + imag_sum * imag_sum);
                    if (magnitude > max_bin_magnitude) {
                        max_bin_magnitude = magnitude;
                    }
                }
                
                // Accumulate the strongest magnitude found for this pitch in this octave
                total_pitch_magnitude += max_bin_magnitude;
            }
            
            // Output layout: [C_L, C_R, C#_L, C#_R, D_L, D_R...]
            d_chroma[pitch_class * num_channels + ch] = total_pitch_magnitude / 6.0f; // Average across the 6 octaves
        }
    }
}

// Host Wrapper
bool compute_chromagram(const int16_t* audio_data, uint32_t num_frames, uint16_t num_channels, float* out_chroma, uint32_t* out_num_pitches) {
    if (!audio_data || !out_chroma || num_frames == 0 || num_channels < 1 || num_channels > 2) return false;

    if (out_num_pitches) {
        *out_num_pitches = 12; // 12 chromatic pitches
    }

    if (current_buffer_frames_chroma != num_frames || current_channels_chroma != num_channels) {
        if (d_audio_buffer_chroma) cudaFree(d_audio_buffer_chroma);
        if (d_chroma_buffer) cudaFree(d_chroma_buffer);
        
        cudaMalloc((void**)&d_audio_buffer_chroma, num_frames * num_channels * sizeof(int16_t));
        cudaMalloc((void**)&d_chroma_buffer, 12 * num_channels * sizeof(float));
        
        current_buffer_frames_chroma = num_frames;
        current_channels_chroma = num_channels;
    }

    cudaMemcpy(d_audio_buffer_chroma, audio_data, num_frames * num_channels * sizeof(int16_t), cudaMemcpyHostToDevice);

    int blockSize = 12; // One block of 12 threads for the 12 pitches is enough
    int numBlocks = 1;
    size_t shared_mem_size = num_frames * num_channels * sizeof(float);
    compute_stereo_chromagram<<<numBlocks, blockSize, shared_mem_size>>>(d_audio_buffer_chroma, d_chroma_buffer, num_frames, num_channels, 44100.0f);
    
    cudaDeviceSynchronize();

    cudaMemcpy(out_chroma, d_chroma_buffer, 12 * num_channels * sizeof(float), cudaMemcpyDeviceToHost);

    return true;
}

void cleanup_gpu_dsp() {
    if (d_audio_buffer_bands) { cudaFree(d_audio_buffer_bands); d_audio_buffer_bands = NULL; }
    if (d_magnitude_buffer_bands) { cudaFree(d_magnitude_buffer_bands); d_magnitude_buffer_bands = NULL; }
    current_buffer_frames_bands = 0;
    current_channels_bands = 0;
    
    if (d_audio_buffer_chroma) { cudaFree(d_audio_buffer_chroma); d_audio_buffer_chroma = NULL; }
    if (d_chroma_buffer) { cudaFree(d_chroma_buffer); d_chroma_buffer = NULL; }
    current_buffer_frames_chroma = 0;
    current_channels_chroma = 0;
}
