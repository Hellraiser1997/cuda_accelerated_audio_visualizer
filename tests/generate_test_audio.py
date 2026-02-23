import wave
import struct
import math
import os

SAMPLE_RATE = 44100
DURATION = 5.0 # seconds

def generate_wav(filename, audio_data, channels=1):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    with wave.open(filepath, 'w') as wav_file:
        wav_file.setnchannels(channels)
        wav_file.setsampwidth(2) # 16-bit
        wav_file.setframerate(SAMPLE_RATE)
        
        for sample in audio_data:
            # Clamp to 16-bit int range
            sample = max(-32768, min(32767, int(sample)))
            data = struct.pack('<h', sample)
            wav_file.writeframesraw(data)
    print(f"Generated: {filename}")

def sine_wave(freq, amplitude=0.8):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        val = amplitude * math.sin(2 * math.pi * freq * t)
        data.append(val * 32767)
    return data

def mixed_wave(freqs, amplitudes):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        val = 0
        for f, a in zip(freqs, amplitudes):
            val += a * math.sin(2 * math.pi * f * t)
        data.append(val * 32767)
    return data

def frequency_sweep(start_freq, end_freq, amplitude=0.8):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        # Exponential sweep formula for more natural progression
        # current_freq = start_freq * (end_freq / start_freq) ** (t / DURATION)
        # Linear sweep is easier to visualize moving left to right
        current_freq = start_freq + (end_freq - start_freq) * (t / DURATION)
        
        # Phase needs to be integrated over time, not just freq * t
        # Integral of (start + rt) = start*t + 0.5*r*t^2
        r = (end_freq - start_freq) / DURATION
        phase = start_freq * t + 0.5 * r * t * t
        
        val = amplitude * math.sin(2 * math.pi * phase)
        data.append(val * 32767)
    return data

def amplitude_sweep(freq):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        amplitude = t / DURATION # Goes from 0.0 to 1.0
        val = amplitude * math.sin(2 * math.pi * freq * t)
        data.append(val * 32767)
    return data

def alternating_frequencies(freq1, freq2, interval=0.5, amplitude=0.8):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        if int(t / interval) % 2 == 0:
            val = amplitude * math.sin(2 * math.pi * freq1 * t)
        else:
            val = amplitude * math.sin(2 * math.pi * freq2 * t)
        data.append(val * 32767)
    return data

def stereo_split_test(freq_left, freq_right, amplitude=0.8):
    data = []
    for i in range(int(SAMPLE_RATE * DURATION)):
        t = i / SAMPLE_RATE
        val_l = amplitude * math.sin(2 * math.pi * freq_left * t)
        val_r = amplitude * math.sin(2 * math.pi * freq_right * t)
        data.append(val_l * 32767)
        data.append(val_r * 32767)
    return data

def run_all():
    # 1. Constant Frequency/Amplitude (Mono)
    generate_wav("01_constant_bass_100hz.wav", sine_wave(100))
    generate_wav("02_constant_mid_1000hz.wav", sine_wave(1000))
    generate_wav("03_constant_high_10000hz.wav", sine_wave(10000))
    
    # 2. Amplitude Sweep
    generate_wav("04_amplitude_sweep_440hz.wav", amplitude_sweep(440))
    
    # 3. Frequency Sweep (Chirp)
    generate_wav("05_frequency_sweep_low_to_high.wav", frequency_sweep(20, 15000))
    
    # 4. Mixed Signals
    generate_wav("06_mixed_signal_chord.wav", mixed_wave([440, 554, 659], [0.33, 0.33, 0.33])) # A Major Chord
    
    # 5. Alternating Test
    generate_wav("07_alternating_bass_treble.wav", alternating_frequencies(100, 8000))

    # 6. Stereo Test (Left Ear Bass, Right Ear Treble)
    generate_wav("08_stereo_split_100L_5000R.wav", stereo_split_test(100, 5000), channels=2)

if __name__ == "__main__":
    run_all()
