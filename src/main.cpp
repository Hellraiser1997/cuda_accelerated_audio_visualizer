#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_ttf.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

extern "C" {
    #include "media_reader.h"
    #include "audio_player.h"
}
#include "gpu_dsp.h"

#define MAX_DISPLAY_BANDS 512

// Helper to draw text via OpenGL Textures for Stats For Nerds HUD
void render_text_ortho(const char* text, float x, float y, TTF_Font* font, SDL_Color color, int win_w, int win_h) {
    if (!text || !font) return;
    SDL_Surface* surface = TTF_RenderText_Blended(font, text, color);
    if (!surface) return;
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    int mode = GL_RGB;
    if(surface->format->BytesPerPixel == 4) {
        mode = (surface->format->Rmask == 0x000000ff) ? GL_RGBA : GL_BGRA;
    }
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, surface->w, surface->h, 0, mode, GL_UNSIGNED_BYTE, surface->pixels);
    
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    gluOrtho2D(0, win_w, win_h, 0);
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();
    
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(x, y);
    glTexCoord2f(1.0f, 0.0f); glVertex2f(x + surface->w, y);
    glTexCoord2f(1.0f, 1.0f); glVertex2f(x + surface->w, y + surface->h);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(x, y + surface->h);
    glEnd();
    
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
    glDeleteTextures(1, &texture);
    SDL_FreeSurface(surface);
}

enum VisMode {
    MODE_BARS = 0,
    MODE_RADIAL = 1,
    MODE_CHROMAGRAM = 2,
    MODE_SPECTROGRAM = 3,
    MODE_OSCILLOSCOPE = 4,
    MODE_PHASE_SPACE = 5,
    MODE_3D_MOUNTAIN = 6,
    MODE_PARTICLE_CLOUD = 7,  // New 3D Spatial Mode
    MODE_3D_LISSAJOUS = 8,    // New 3D Spatial Mode
    MODE_ULTIMATE_GALAXY = 9, // Ultimate Visualizer
    MODE_COUNT = 10           // Determines wrap limit
};

// Struct for Mode 2: Stereo Positional Ripples
typedef struct {
    float pan_x;      // Stereo position (0.0 Left ... 1.0 Right)
    float radius;     // Expanding radius of the ripple
    float intensity;  // How hard the wave hit
    float life;       // 1.0 down to 0.0 fade out
    bool active;
} Ripple;

// Physics struct for Volumetric Audio Clouds
struct Particle {
    float x, y, z;
    float vx, vy, vz;
    float life; // Used to respawn dead particles dynamically
};


int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s <input.wav>\n", argv[0]);
        return 1;
    }

    const char* filepath = argv[1];
    
    // 1. Load Audio Media via FFmpeg
    printf("Loading Audio Media file: %s\n", filepath);
    AudioData audio_data;
    if (!load_media(filepath, &audio_data)) {
        return 1;
    }
    
    printf("Loaded. Frames: %u. Channels: %u. Sample Rate: %u\n", audio_data.num_samples, audio_data.header.num_channels, audio_data.header.sample_rate);
    
    // 2. Initialize GPU
    if (!init_gpu_dsp()) {
        fprintf(stderr, "Failed to initialize CUDA GPU DSP.\n");
        free_media(&audio_data);
        return 1;
    }
    
    // 3. Initialize Audio Player
    AudioPlayer player;
    if (!init_audio_player(&player, &audio_data)) {
        fprintf(stderr, "Failed to initialize audio player.\n");
        free_media(&audio_data);
        cleanup_gpu_dsp();
        return 1;
    }

    // 4. Initialize SDL for Graphics
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        cleanup_audio_player(&player);
        cleanup_gpu_dsp();
        free_media(&audio_data);
        return 1;
    }

    int window_width = 1000;
    int window_height = 800;
    
    // Request an OpenGL 2.1 Context
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    
    SDL_Window* window = SDL_CreateWindow("CUDA Music Player (OpenGL 3D) - High Res Interactive GPU Equalizer",
                                          SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                          window_width, window_height,
                                          SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if (!window) {
        fprintf(stderr, "Window could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_Quit();
        cleanup_audio_player(&player);
        cleanup_gpu_dsp();
        free_media(&audio_data);
        return 1;
    }

    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (!gl_context) {
        fprintf(stderr, "OpenGL context could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        cleanup_audio_player(&player);
        cleanup_gpu_dsp();
        free_media(&audio_data);
        return 1;
    }
    
    if (TTF_Init() == -1) {
        fprintf(stderr, "TTF_Init failed: %s\n", TTF_GetError());
        return 1;
    }
    TTF_Font* ui_font = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24);
    if (!ui_font) {
        fprintf(stderr, "Failed to load fonts, text overlay will not work.\n");
    }
    
    // Enable VSync
    SDL_GL_SetSwapInterval(1);
    
    // Basic OpenGL Init
    glClearColor(20.0f/255.0f, 20.0f/255.0f, 25.0f/255.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // 5. Start Playback
    printf("Starting playback...\n");
    start_playback(&player);

    uint32_t num_frames_per_chunk = 8192; // Ultra High resolution sampling window
    uint16_t num_channels = audio_data.header.num_channels;
    
    // Arrays for the GPU to return data into. 
    float* magnitudes = (float*)malloc(MAX_DISPLAY_BANDS * 2 * sizeof(float)); 
    uint32_t display_bands_returned;
    
    float* chromagram = (float*)malloc(12 * 2 * sizeof(float));
    uint32_t chroma_pitches;
    
    // Setup for 3D Mountain Range Buffer (100 frames deep x MAX_DISPLAY_BANDS wide x 2 channels)
    int mountain_depth = 100;
    float** mountain_history = (float**)malloc(mountain_depth * sizeof(float*));
    for (int i=0; i<mountain_depth; ++i) {
        mountain_history[i] = (float*)calloc(MAX_DISPLAY_BANDS * 2, sizeof(float));
    }
    
    // Setup for Volumetric Particle Cloud
    int num_particles = 10000;
    Particle* particles = (Particle*)calloc(num_particles, sizeof(Particle));
    for (int i = 0; i < num_particles; ++i) {
        // Spawn randomly in a large bounding box
        particles[i].x = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].y = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].z = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].life = (rand() % 100) / 100.0f;
    }

    // Setup for Ultimate Procedural Galaxy
    int num_galaxy_particles = 20000;
    Particle* galaxy_particles = (Particle*)calloc(num_galaxy_particles, sizeof(Particle));
    for (int i = 0; i < num_galaxy_particles; ++i) {
        float angle = ((float)rand() / RAND_MAX) * 20.0f * M_PI; // Spirals
        float radius = ((float)rand() / RAND_MAX) * 250.0f + (angle * 4.0f);
        galaxy_particles[i].x = cosf(angle) * radius;
        galaxy_particles[i].y = ((float)rand() / RAND_MAX) * 20.0f - 10.0f; // Thickness
        galaxy_particles[i].z = sinf(angle) * radius;
        galaxy_particles[i].vx = angle; // store original angle for animation
        galaxy_particles[i].vy = radius; // store original radius
        galaxy_particles[i].vz = ((float)rand() / RAND_MAX) * M_PI * 2.0f; // random phase
    }
    // Setup for Straight Line Grid (Mode 1 Replacement)
    int num_lines = 600;
    int points_per_line = 1200;
    int total_grid_particles = num_lines * points_per_line;
    Particle* grid_particles = (Particle*)calloc(total_grid_particles, sizeof(Particle));
    
    // Create a 2D plane grid of straight horizontal lines
    float grid_size_x = 900.0f;
    float grid_size_z = 900.0f;
    
    for (int l = 0; l < num_lines; ++l) {
        float z = -grid_size_z/2.0f + ((float)l / (num_lines - 1)) * grid_size_z;
        for (int p = 0; p < points_per_line; ++p) {
            float x = -grid_size_x/2.0f + ((float)p / (points_per_line - 1)) * grid_size_x;
            int idx = l * points_per_line + p;
            
            grid_particles[idx].x = x;
            grid_particles[idx].y = 0.0f; 
            grid_particles[idx].z = z;
            grid_particles[idx].vx = 0.0f;         // Unused
            grid_particles[idx].vy = 0.0f;         // Unused
            grid_particles[idx].vz = 0.0f;         // Spring velocity
            grid_particles[idx].life = 0.0f;       // Unused
        }
    }

    const char* mode_names[] = {
        "CUDA Music Player - Mode 0: Stereo Logarithmic Bar EQ",
        "CUDA Music Player - Mode 1: [3D] Particle Ring Vibrator",
        "CUDA Music Player - Mode 2: Stereo Positional Waves (Digital Pond)",
        "CUDA Music Player - Mode 3: Spectrogram Waterfall",
        "CUDA Music Player - Mode 4: The Oscilloscope (Vectorscope)",
        "CUDA Music Player - Mode 5: Phase Space Attractor",
        "CUDA Music Player - Mode 6: [3D] Cyber-Mountain Terrain",
        "CUDA Music Player - Mode 7: [3D] Volumetric Particle Cloud",
        "CUDA Music Player - Mode 8: [3D] Lissajous Phase Sphere",
        "CUDA Music Player - Mode 9: [3D] Procedural Audio Galaxy (Ultimate)"
    };
    int last_mode_index = -1;

    // Visualization Loop
    bool quit = false;
    SDL_Event e;
    
    int current_mode_index = 0;
    bool is_fullscreen = false;
    bool show_stats = false;
    uint32_t last_time = SDL_GetTicks();
    int frames = 0;
    float current_fps = 0.0f;
    
    GLUquadric* planet_quad = gluNewQuadric();

    while (player.is_playing && !quit) {
        if (current_mode_index != last_mode_index) {
            SDL_SetWindowTitle(window, mode_names[current_mode_index]);
            last_mode_index = current_mode_index;
        }

        // Handle window events
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
                player.is_playing = false; // Stop audio thread
            }
            if (e.type == SDL_WINDOWEVENT) {
                if (e.window.event == SDL_WINDOWEVENT_RESIZED) {
                    window_width = e.window.data1;
                    window_height = e.window.data2;
                    glViewport(0, 0, window_width, window_height);
                }
            }
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_f) {
                    is_fullscreen = !is_fullscreen;
                    SDL_SetWindowFullscreen(window, is_fullscreen ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
                }
                else if (e.key.keysym.sym == SDLK_s) {
                    show_stats = !show_stats;
                }
                // Number Key Shortcuts for Modes 0-8
                else if (e.key.keysym.sym >= SDLK_1 && e.key.keysym.sym <= SDLK_9) {
                    current_mode_index = e.key.keysym.sym - SDLK_1;
                }
                // Unlimited Mode Cycling Navigation
                else if (e.key.keysym.sym == SDLK_RIGHT) {
                    current_mode_index = (current_mode_index + 1) % MODE_COUNT;
                }
                else if (e.key.keysym.sym == SDLK_LEFT) {
                    current_mode_index = (current_mode_index - 1 + MODE_COUNT) % MODE_COUNT;
                }
                else if (e.key.keysym.sym == SDLK_SPACE) {
                    player.is_playing = !player.is_playing;
                    if (player.is_playing) start_playback(&player);
                }
            }
        }

        uint32_t play_index = player.current_sample_index; // Measured in frames
        
        if (play_index + num_frames_per_chunk < audio_data.num_samples) {
            
            // `play_index` is frames, since audio_data is a flat array of int16_t, 
            // the offset in the array is `play_index * num_channels`
            size_t sample_offset = play_index * num_channels;
            
            // Get frequency magnitudes from GPU (Now doing Stereo and Logarithmic Bucketing!)
            process_audio_chunk(&audio_data.audio_data[sample_offset], num_frames_per_chunk, num_channels, magnitudes, &display_bands_returned);
            int display_bands = display_bands_returned;
            
            // Shift mountain history and add new magnitudes
            for (int r = mountain_depth - 1; r > 0; --r) {
                memcpy(mountain_history[r], mountain_history[r-1], 200 * 2 * sizeof(float));
            }
            memcpy(mountain_history[0], magnitudes, 200 * 2 * sizeof(float));
            
            // Audio analysis for global parameters
            float rms = 0;
            float bass_sum = 0;
            float peak_mag = 0;
            int peak_band = 0;
            for(int i = 0; i < display_bands; ++i) {
                float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                rms += mag * mag;
                if(i < 5) bass_sum += mag;
                if(mag > peak_mag) { peak_mag = mag; peak_band = i; }
            }
            rms = sqrtf(rms / display_bands);

            // Wait to clear the target frame entirely unless we are in the spectrogram mode which does its own target management
            if (current_mode_index != MODE_SPECTROGRAM) {
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            }
            if (current_mode_index == MODE_BARS) {
                // =============== BAR VISUALIZER (STEREO SPLIT) ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluOrtho2D(0, window_width, window_height, 0); 
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                int bar_width = window_width / display_bands;
                int padding = 2; // Space between bars
                
                glBegin(GL_QUADS);
                for (int band = 0; band < display_bands; ++band) {
                    
                    // Channel 0 (Left), Channel 1 (Right) if stereo.
                    // If mono, both use index 0.
                    float mag_left = magnitudes[band * num_channels + 0];
                    float mag_right = magnitudes[band * num_channels + (num_channels > 1 ? 1 : 0)];
                    
                    // Scale magnitude (Logarithmic processing means magnitudes are much larger/denser generally)
                    int h_left = (int)(mag_left * 2.0f); 
                    int h_right = (int)(mag_right * 2.0f);
                    
                    // Floor and ceil
                    if (h_left > window_height / 2) h_left = (window_height / 2) - 10;
                    if (h_left < 2) h_left = 2;
                    if (h_right > window_height / 2) h_right = (window_height / 2) - 10;
                    if (h_right < 2) h_right = 2;
                    
                    // Color gradient logic (Blue to Pink)
                    Uint8 r = (band * 255) / display_bands;
                    Uint8 g = 100 + (r / 4);
                    Uint8 b = 255 - ((band * 128) / display_bands);
                    glColor3ub(r, g, b);
                    
                    // Draw Left Channel (Growing UP from the center line)
                    int h_center = window_height / 2;
                    int x = band * bar_width + padding;
                    int bw = bar_width - (padding * 2);
                    
                    glVertex2f(x, h_center);
                    glVertex2f(x + bw, h_center);
                    glVertex2f(x + bw, h_center - h_left);
                    glVertex2f(x, h_center - h_left);
                    
                    // Draw Right Channel (Growing DOWN from the center line)
                    glVertex2f(x, h_center);
                    glVertex2f(x + bw, h_center);
                    glVertex2f(x + bw, h_center + h_right);
                    glVertex2f(x, h_center + h_right);
                }
                glEnd();
            } else if (current_mode_index == MODE_RADIAL) {
                // =============== [3D] PARTICLE RING VIBRATOR (Refined with Physics) ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity();
                // Lower FOV for a more cinematic, deep look
                gluPerspective(50.0f, (float)window_width / (float)window_height, 0.1f, 1500.0f);
                
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                // Audio analysis for Global parameters
                float bass_sum = 0;
                float current_frame_max = 0.0f;
                
                for(int i = 0; i < (int)(display_bands * (5.0f/80.0f)); ++i) {
                    float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                    bass_sum += mag;
                    if (mag > current_frame_max) current_frame_max = mag;
                }
                float mid_sum = 0;
                for(int i=(int)(display_bands * 0.25f); i<(int)(display_bands * 0.5f); ++i) {
                     float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                     mid_sum += mag;
                     if (mag > current_frame_max) current_frame_max = mag;
                }
                
                // Check all bands for AGC
                for (int i=0; i<display_bands; ++i) {
                    float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                    if (mag > current_frame_max) current_frame_max = mag;
                }
                
                // AGC Logic (Automatic Gain Control)
                static float rolling_peak = 1.0f;
                if (current_frame_max > rolling_peak) {
                    rolling_peak = current_frame_max; // Instant rise for punch
                } else {
                    rolling_peak *= 0.995f; // Smooth decay back down
                    if (rolling_peak < 1.0f) rolling_peak = 1.0f; // Minimum floor to prevent amplifying silence
                }
                
                // Calculate dynamic scale factor
                float target_peak = 2.0f; // Aesthetically pleasing target max magnitude
                float dynamic_scale_factor = target_peak / rolling_peak;
                
                // Scale global trackers
                bass_sum *= dynamic_scale_factor;
                mid_sum *= dynamic_scale_factor;
                
                // Position Camera lower and looking across the horizon for depth
                static float ring_cam_rot = 0;
                ring_cam_rot += 0.10f; // Smooth, slow pan
                
                // Keep the camera stable, only subtle bumping, moved closer to the grid for depth.
                float cam_dist = 180.0f + (bass_sum * 1.5f); 
                
                gluLookAt(sinf(ring_cam_rot * 0.01f) * cam_dist, 35.0f - (bass_sum * 1.0f), cosf(ring_cam_rot * 0.01f) * cam_dist, 
                          0.0f, -10.0f, 0.0f, 
                          0.0f, 1.0f, 0.0f);

                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Additive glowing blending
                
                // Massive Geometry Bloom setup
                glEnable(GL_POINT_SMOOTH); 
                glPointSize(1.0f); // Minimum pixel size for maximum point density

                // Time factor for standing waves
                float t = SDL_GetTicks() * 0.003f;
                
                // --- GRID WITH FLUID PHYSICS (Lines) ---
                glBegin(GL_POINTS);
                for (int l = 0; l < num_lines; ++l) {
                    for (int p = 0; p < points_per_line; ++p) {
                        int i = l * points_per_line + p;
                        Particle* pt = &grid_particles[i];
                        
                        // Radial Audio Displacement:
                        float radius = sqrtf(pt->x * pt->x + pt->z * pt->z);
                        
                        // Map the full 512 visually active bands across the expansive 900x900 grid
                        float max_radius_for_freq = 400.0f;
                        int band = (int)((radius / max_radius_for_freq) * display_bands);
                        if (band >= display_bands) band = display_bands - 1;
                        if (band < 0) band = 0;
                        
                        // mapped history: the wave moves outwards from the center over time
                        int history_idx = (int)((radius / max_radius_for_freq) * mountain_depth);
                        if (history_idx < 0) history_idx = 0;
                        if (history_idx >= mountain_depth) history_idx = mountain_depth - 1;
                        
                        float mag = (mountain_history[history_idx][band * num_channels + 0] + 
                                     mountain_history[history_idx][band * num_channels + (num_channels>1?1:0)]) / 2.0f;
                        
                        // Apply AGC Dynamic Scale
                        mag *= dynamic_scale_factor;
                                     
                        // Deep funnel valley shape
                        float norm_r = radius / max_radius_for_freq;
                        if (norm_r > 1.0f) norm_r = 1.0f;
                        float bowl_shape = (powf(norm_r, 3.0f) * 120.0f) - 60.0f; 

                        // Smooth cone displacement based on bass
                        float cone_displacement = 0.0f;
                        if (radius < 100.0f) {
                           cone_displacement = bass_sum * (100.0f - radius) * 0.4f; 
                        }
                        
                        // Controlled Ripples
                        float ripple_height = mag * 50.0f; 
                        
                        // Subtle standing wave 
                        float standing_wave = sinf(radius * 0.1f - t) * (norm_r) * mag * 8.0f;

                        float target_y = bowl_shape + cone_displacement + ripple_height + standing_wave;

                        // SPRING PHYSICS
                        float stiffness = 0.2f; 
                        float damping = 0.80f;
                        
                        pt->vz += (target_y - pt->y) * stiffness; 
                        pt->vz *= damping;
                        pt->y += pt->vz;

                        // Clean, Silver/White aesthetic
                        float intensity = mag * 1.5f + 0.15f;
                        if (intensity > 1.0f) intensity = 1.0f;
                        
                        Uint8 col_r = (Uint8)(255 * intensity);
                        Uint8 col_g = (Uint8)(255 * intensity);
                        Uint8 col_b = (Uint8)(255 * intensity);

                        if (intensity < 0.25f) {
                             col_r *= 0.7f;
                             col_g *= 0.8f;
                        }

                        // Strict Red transitions: ONLY for the highest peaks
                        if (pt->y > 20.0f && radius > 150.0f) {
                            float peak_factor = (pt->y - 20.0f) / 40.0f; 
                            if (peak_factor > 1.0f) peak_factor = 1.0f;
                            
                            col_r = 255;
                            col_g = (Uint8)(255 * (1.0f - peak_factor)); 
                            col_b = (Uint8)(255 * (1.0f - peak_factor));
                        }
                        
                        // Simulated BLOOM via very dense stacking with low alpha
                        Uint8 alpha = (Uint8)(intensity * 22.0f); // Exceptionally low base alpha (around ~10/255)
                        if (alpha > 255) alpha = 255;
                        
                        glColor4ub(col_r, col_g, col_b, alpha);
                        glVertex3f(pt->x, pt->y, pt->z);
                    }
                }
                glEnd();
                
                glDisable(GL_BLEND);
                glDisable(GL_POINT_SMOOTH);
                glDisable(GL_BLEND);
                glDisable(GL_POINT_SMOOTH);
            } else if (current_mode_index == MODE_CHROMAGRAM) {
                // =============== MODE 2: STEREO POSITIONAL WAVES (Digital Pond) ===============
                // Replaced Chromagram with the Positional Wave matrix
                glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluOrtho2D(0, window_width, window_height, 0);
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                // Static array of ripples to manage the physics of the pond
                #define MAX_RIPPLES 100
                static Ripple ripples[MAX_RIPPLES] = {0};
                
                // 1. Spawning Logic: Analyze audio and drop new ripples
                int num_macro_bands = 16;
                int bands_per_macro = display_bands / num_macro_bands;
                
                for (int m = 0; m < num_macro_bands; ++m) {
                    float sum_L = 0;
                    float sum_R = 0;
                    
                    // Sum energy in this macro band to prevent chaotic overlapping
                    for (int b = 0; b < bands_per_macro; ++b) {
                        int actual_band = m * bands_per_macro + b;
                        // Magnitudes are interleaved [L, R, L, R...]
                        sum_L += magnitudes[actual_band * num_channels + 0];
                        if (num_channels > 1) {
                            sum_R += magnitudes[actual_band * num_channels + 1];
                        } else {
                            sum_R += magnitudes[actual_band * num_channels + 0]; // fallback to mono
                        }
                    }
                    
                    float total_energy = sum_L + sum_R;
                    
                    // If there's a strong beat in this macro band, drop a ripple!
                    if (total_energy > 0.8f) { 
                        float pan_x = (total_energy > 0.001f) ? (sum_R / total_energy) : 0.5f;
                        
                        // Find a dead ripple slot to recycle
                        for (int i = 0; i < MAX_RIPPLES; ++i) {
                            if (!ripples[i].active) {
                                ripples[i].active = true;
                                ripples[i].pan_x = pan_x;
                                ripples[i].radius = 1.0f;
                                ripples[i].life = 1.0f;
                                ripples[i].intensity = total_energy;
                                break;
                            }
                        }
                    }
                }
                
                // 2. Physics & Rendering Logic
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Glowing additive blend
                glEnable(GL_LINE_SMOOTH);
                glLineWidth(2.0f);
                
                for (int i = 0; i < MAX_RIPPLES; ++i) {
                    if (ripples[i].active) {
                        // Expand the circle based on how hard it hit
                        ripples[i].radius += ripples[i].intensity * 5.0f + 2.0f;
                        ripples[i].life -= 0.02f; // Fade out gradually
                        
                        if (ripples[i].life <= 0.0f || ripples[i].radius > window_width) {
                            ripples[i].active = false;
                            continue;
                        }
                        
                        // Map the physical 0.0-1.0 stereo pan to the X-axis of the screen
                        float center_x = ripples[i].pan_x * window_width;
                        float center_y = window_height / 2.0f; // Drops hit the center line (Y-axis)
                        
                        // Calculate color based on Life and Intensity
                        Uint8 alpha = (Uint8)(ripples[i].life * 200.0f);
                        if (alpha > 255) alpha = 255;
                        
                        // Let's color code it. Left = Cyan, Center = White, Right = Magenta
                        Uint8 r = (Uint8)(ripples[i].pan_x * 255);
                        Uint8 g = (Uint8)((1.0f - fabs(ripples[i].pan_x - 0.5f) * 2.0f) * 200) + 55;
                        Uint8 b = (Uint8)((1.0f - ripples[i].pan_x) * 255);
                        
                        glColor4ub(r, g, b, alpha);
                        
                        // Draw the expanding circular ripple (Using a Line Loop)
                        glBegin(GL_LINE_LOOP);
                        int num_segments = 60;
                        for (int seg = 0; seg < num_segments; ++seg) {
                            float theta = 2.0f * 3.1415926f * (float)seg / (float)num_segments;
                            float x = center_x + ripples[i].radius * cosf(theta);
                            float y = center_y + ripples[i].radius * sinf(theta);
                            glVertex2f(x, y);
                        }
                        glEnd();
                    }
                }
                
                glDisable(GL_BLEND);
                glDisable(GL_LINE_SMOOTH);
                
            } else if (current_mode_index == MODE_SPECTROGRAM) {
                // =============== SPECTROGRAM (Orthographic Waterfall) ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluOrtho2D(0, window_width, window_height, 0); 
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                int col_width = window_width / display_bands;
                int row_height = window_height / mountain_depth;
                if (row_height < 1) row_height = 1;
                
                glBegin(GL_QUADS);
                for (int row = 0; row < mountain_depth; ++row) {
                    for (int band = 0; band < display_bands; ++band) {
                        float mag_left = mountain_history[row][band * num_channels + 0];
                        float mag_right = mountain_history[row][band * num_channels + (num_channels > 1 ? 1 : 0)];
                        float avg_mag = (mag_left + mag_right) / 2.0f;
                        
                        int intensity = (int)(avg_mag * 1200);
                        if (intensity > 255) intensity = 255;
                        
                        glColor3ub(intensity, intensity / 2, intensity / 4);
                        
                        int x = band * col_width;
                        int y = row * row_height;
                        glVertex2f(x, y + row_height);
                        glVertex2f(x + col_width, y + row_height);
                        glVertex2f(x + col_width, y);
                        glVertex2f(x, y);
                    }
                }
                glEnd();

            } else if (current_mode_index == MODE_OSCILLOSCOPE) {
                // =============== OSCILLOSCOPE (Lissajous) ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluOrtho2D(0, window_width, window_height, 0); 
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                glColor3ub(50, 255, 100); // Glowing Green
                int center_x = window_width / 2;
                int center_y = window_height / 2;
                
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < num_frames_per_chunk; ++i) {
                     int16_t l_curr = audio_data.audio_data[sample_offset + i * num_channels + 0];
                     int16_t r_curr = audio_data.audio_data[sample_offset + i * num_channels + (num_channels > 1 ? 1 : 0)];
                     
                     int x2 = center_x + (int)((r_curr / 32768.0f) * (window_width / 2) * 0.8f);
                     int y2 = center_y - (int)((l_curr / 32768.0f) * (window_height / 2) * 0.8f);
                     
                     glVertex2f(x2, y2);
                }
                glEnd();
                
            } else if (current_mode_index == MODE_PHASE_SPACE) {
                // =============== PHASE SPACE ATTRACTOR ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluOrtho2D(0, window_width, window_height, 0); 
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();
                
                glColor3ub(255, 100, 255); // Pink
                int center_x = window_width / 2;
                int center_y = window_height / 2;
                int delay_samples = 20; // Tau
                
                glBegin(GL_LINE_STRIP);
                for (int i = delay_samples; i < num_frames_per_chunk; ++i) {
                     float p_curr = (audio_data.audio_data[sample_offset + i * num_channels + 0] + audio_data.audio_data[sample_offset + i * num_channels + (num_channels > 1 ? 1 : 0)]) / 2.0f;
                     float p_delay_curr = (audio_data.audio_data[sample_offset + (i-delay_samples) * num_channels + 0] + audio_data.audio_data[sample_offset + (i-delay_samples) * num_channels + (num_channels > 1 ? 1 : 0)]) / 2.0f;
                     
                     int x2 = center_x + (int)((p_curr / 32768.0f) * (window_width / 2) * 0.8f);
                     int y2 = center_y + (int)((p_delay_curr / 32768.0f) * (window_height / 2) * 0.8f);
                     
                     glVertex2f(x2, y2);
                }
                glEnd();
            } else if (current_mode_index == MODE_3D_MOUNTAIN) {
                // =============== 3D MOUNTAIN TERRAIN ===============
                // Set true 3D perspective camera matrix
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f, (float)window_width / (float)window_height, 0.1f, 1000.0f);
                
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();
                
                // Position Camera slightly above and looking down
                gluLookAt(0.0f, 40.0f, -60.0f,    // Camera position
                          0.0f, -10.0f,  30.0f,    // Target position
                          0.0f, 1.0f,   0.0f);     // Up vector
                          
                // Spin it very slowly automatically 
                static float rot = 0;
                rot += 0.2f;
                glRotatef(rot, 0, 1, 0);

                // Draw the mountain terrain wireframe using the history buffer
                glBegin(GL_LINES);
                float width_spread = 160.0f; // Total width of the grid 
                float depth_spacing = 1.2f;  // Space between rows
                
                for(int row = 0; row < mountain_depth - 1; row++) {
                    for(int col = 0; col < display_bands - 1; col++) {
                        
                        // Grab magnitude averages from buffer
                        float m1 = (mountain_history[row][col*num_channels+0] + mountain_history[row][col*num_channels+(num_channels>1?1:0)])/2.0f;
                        float m2 = (mountain_history[row][(col+1)*num_channels+0] + mountain_history[row][(col+1)*num_channels+(num_channels>1?1:0)])/2.0f;
                        float m3 = (mountain_history[row+1][col*num_channels+0] + mountain_history[row+1][col*num_channels+(num_channels>1?1:0)])/2.0f;
                        
                        float x1 = (col / (float)display_bands) * width_spread - (width_spread/2.0f);
                        float x2 = ((col+1) / (float)display_bands) * width_spread - (width_spread/2.0f);
                        
                        float z1 = row * depth_spacing;
                        float z2 = (row+1) * depth_spacing;
                        
                        float y1 = m1 * 400.0f; // Scale height mapping
                        float y2 = m2 * 400.0f;
                        float y3 = m3 * 400.0f;
                        
                        // Synthwave Neon Cyan shader
                        glColor3ub(0, 200, 255);
                        
                        // Draw Vertical connection
                        glVertex3f(x1, y1, z1);
                        glVertex3f(x1, y3, z2);
                    }
                }
                glEnd();
            } else if (current_mode_index == MODE_PARTICLE_CLOUD) {
                // =============== VOLUMETRIC PARTICLE CLOUD ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity();
                gluPerspective(60.0f, (float)window_width / (float)window_height, 0.1f, 1000.0f);
                glMatrixMode(GL_MODELVIEW);  glLoadIdentity();
                
                gluLookAt(0.0f, 0.0f, -150.0f,    // Camera position
                          0.0f, 0.0f,  0.0f,      // Look at center
                          0.0f, 1.0f,  0.0f);     // Up vector
                
                // Spin slowly to see the 3D depth of the cloud
                static float rot_p = 0;
                rot_p += 0.5f;
                glRotatef(rot_p, 0, 1, 0);
                glRotatef(rot_p * 0.3f, 1, 0, 0); // Slight wobble pitch

                // Extract audio dynamics to influence physics
                float bass_sum = 0;   // Low frequencies (gravity well / red)
                float treble_sum = 0; // High frequencies (explosive outward force / blue)
                
                for(int i = 0; i < 5; ++i) bass_sum += magnitudes[i * num_channels];
                for(int i = display_bands - 10; i < display_bands; ++i) treble_sum += magnitudes[i * num_channels];
                bass_sum /= 5.0f;
                treble_sum /= 10.0f;

                glBegin(GL_POINTS);
                for (int i = 0; i < num_particles; ++i) {
                    Particle* p = &particles[i];

                    // 1. Gravity Well (Pulling toward 0,0,0 proportional to Bass)
                    float dist_sq = (p->x * p->x) + (p->y * p->y) + (p->z * p->z);
                    if (dist_sq < 1.0f) dist_sq = 1.0f; // Avoid singularity explosion
                    
                    float pull_force = (bass_sum * 50.0f) / dist_sq; 
                    p->vx -= (p->x * pull_force);
                    p->vy -= (p->y * pull_force);
                    p->vz -= (p->z * pull_force);

                    // 2. High Frequency Explosion (Pushing outward proportional to Treble)
                    float push_force = treble_sum * 1.5f; 
                    p->vx += (p->x / sqrtf(dist_sq)) * push_force;
                    p->vy += (p->y / sqrtf(dist_sq)) * push_force;
                    p->vz += (p->z / sqrtf(dist_sq)) * push_force;

                    // 3. Simple Friction / Damping
                    p->vx *= 0.95f; p->vy *= 0.95f; p->vz *= 0.95f;

                    // Apply Velocity
                    p->x += p->vx; p->y += p->vy; p->z += p->vz;
                    
                    // Manage Lifespan / Respawn
                    p->life -= 0.01f;
                    if (p->life <= 0 || dist_sq > 60000.0f) {
                        p->x = ((rand() % 4000) - 2000) / 10.0f;
                        p->y = ((rand() % 4000) - 2000) / 10.0f;
                        p->z = ((rand() % 4000) - 2000) / 10.0f;
                        p->vx = p->vy = p->vz = 0.0f;
                        p->life = 1.0f;
                    }

                    // Color mapped dynamically: Red for fast inner particles (bass), Blue for slow outer particles (treble)
                    float speed_color = sqrtf(p->vx*p->vx + p->vy*p->vy + p->vz*p->vz) * 10.0f;
                    if (speed_color > 255) speed_color = 255;
                    
                    glColor3ub((uint8_t)speed_color, (uint8_t)(speed_color/2.0f), (uint8_t)(255 - speed_color));
                    glVertex3f(p->x, p->y, p->z);
                }
                glEnd();

            } else if (current_mode_index == MODE_3D_LISSAJOUS) {
                // =============== 3D LISSAJOUS PHASE SPHERE ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity();
                gluPerspective(60.0f, (float)window_width / (float)window_height, 0.1f, 1000.0f);
                glMatrixMode(GL_MODELVIEW);  glLoadIdentity();
                
                gluLookAt(0.0f, 0.0f, -60.0f,    // Push camera back to see the whole sphere
                          0.0f, 0.0f,   0.0f,    // Target center
                          0.0f, 1.0f,   0.0f);   // Up
                          
                // Spin sphere continuously so we can trace its 3D depth geometry
                static float rot_l = 0;
                rot_l += 0.8f;
                glRotatef(rot_l, 1.0f, 1.0f, 0.0f); 

                int delay_samples = 45; // Approx 1 millisecond at 44.1kHz
                
                glBegin(GL_LINE_STRIP);
                for (int i = delay_samples; i < num_frames_per_chunk; ++i) {
                     int16_t l_curr = audio_data.audio_data[sample_offset + i * num_channels + 0];
                     int16_t r_curr = audio_data.audio_data[sample_offset + i * num_channels + (num_channels > 1 ? 1 : 0)];
                     
                     // Z is a delayed Mono phase composite
                     float p_delay_curr = (audio_data.audio_data[sample_offset + (i-delay_samples) * num_channels + 0] + 
                                           audio_data.audio_data[sample_offset + (i-delay_samples) * num_channels + (num_channels > 1 ? 1 : 0)]) / 2.0f;
                     
                     // Extrapolate 16-bit PCM signals directly into physical 3D vertices
                     float x = (r_curr / 32768.0f) * 40.0f; // Right -> X
                     float y = (l_curr / 32768.0f) * 40.0f; // Left  -> Y
                     float z = (p_delay_curr / 32768.0f) * 40.0f; // Delayed Mono -> Z
                     
                     // Color the knot with a vibrant neon geometry mapping that highlights coordinate density
                     float geom_dist = sqrtf(x*x + y*y + z*z) / 40.0f; 
                     int rgb_int = (int)(geom_dist * 500);
                     if (rgb_int > 255) rgb_int = 255;
                     
                     glColor3ub(50 + rgb_int, 255 - rgb_int, 255);
                     glVertex3f(x, y, z);
                }
                glEnd();
            } else if (current_mode_index == MODE_ULTIMATE_GALAXY) {
                // =============== QUANTUM PARTICLE WAVE MATRIX (ULTIMATE) ===============
                glMatrixMode(GL_PROJECTION); glLoadIdentity();
                gluPerspective(75.0f, (float)window_width / (float)window_height, 0.1f, 1000.0f);
                glMatrixMode(GL_MODELVIEW);  glLoadIdentity();
                
                // Audio analysis for Global parameters
                float rms = 0;
                float bass_sum = 0;
                float peak_mag = 0;
                int peak_band = 0;
                
                for(int i = 0; i < display_bands; ++i) {
                    float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                    rms += mag * mag;
                    if(i < 5) bass_sum += mag;
                    if(mag > peak_mag) { peak_mag = mag; peak_band = i; }
                }
                rms = sqrtf(rms / display_bands);

                // Push camera into the wave matrix and slowly pan/rotate
                static float cam_rot = 0;
                cam_rot += 0.1f + (rms * 0.5f);
                gluLookAt(0.0f, 40.0f, -80.0f,         // Camera moving forward into Z
                          0.0f, -20.0f, 50.0f,         // Looking forward
                          0.0f, 1.0f, 0.0f);           // Up vector
                glRotatef(cam_rot, 0, 1, 0);

                // Rendering 20,000 dense particles as a topographical flow
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Additive glowing blending
                glEnable(GL_POINT_SMOOTH); 
                glPointSize(2.0f); // Thick glowing stars
                
                glBegin(GL_POINTS);
                
                float grid_width = 300.0f;
                float depth_spacing = 3.0f;
                float hue_base = (float)peak_band / display_bands;
                
                // Draw from the back (mountain_depth) to the front
                for(int row = mountain_depth - 1; row >= 0; row--) {
                    for(int col = 0; col < display_bands; col++) {
                        
                        // Extract historical audio magnitude for this exact point in time/frequency
                        float mag = (mountain_history[row][col*num_channels+0] + mountain_history[row][col*num_channels+(num_channels>1?1:0)]) / 2.0f;
                        
                        // X axis is frequency spread
                        float x = (col / (float)display_bands) * grid_width - (grid_width/2.0f);
                        // Z axis is time history
                        float z = row * depth_spacing - (mountain_depth * depth_spacing / 2.0f);
                        // Y axis is geometric height
                        float y = mag * 150.0f; 

                        // Topographical Wave Twisting (The Elegance Math)
                        // Make the edges of the grid curl upward and twist into a ribbon cylinder in deep space
                        float dist_from_center_x = fabsf(x) / (grid_width/2);
                        float curl = powf(dist_from_center_x, 3.0f) * 60.0f; // Edges curl up
                        y += curl; 
                        
                        // Twist the entire grid along the Z-axis
                        float twist_angle = (z * 0.01f) + (SDL_GetTicks() * 0.001f);
                        float tx = x * cosf(twist_angle) - y * sinf(twist_angle);
                        float ty = x * sinf(twist_angle) + y * cosf(twist_angle);
                        x = tx; y = ty;

                        // Color matches Reference Image (Intense Glowing Blues/Greens/Whites)
                        Uint8 r = (Uint8)((sinf((hue_base + (float)col/display_bands) * M_PI * 2.0f) * 0.5f + 0.5f) * 255);
                        Uint8 g = (Uint8)((sinf((hue_base + (float)col/display_bands + 0.33f) * M_PI * 2.0f) * 0.5f + 0.5f) * 255);
                        Uint8 b = (Uint8)((sinf((hue_base + (float)col/display_bands + 0.66f) * M_PI * 2.0f) * 0.5f + 0.5f) * 255);
                        
                        // Alpha/Brightness bound to height (silent valleys are almost invisible, peaks are blindingly bright)
                        float brightness = mag * 3.0f + 0.1f; // base glow
                        if (brightness > 1.0f) brightness = 1.0f;
                        
                        // Bass Flash
                        if (col < 5 && mag > 1.5f) { r=255; g=255; b=255; brightness = 1.0f; }

                        glColor4f(r/255.0f, g/255.0f, b/255.0f, brightness);
                        glVertex3f(x, y, z);
                    }
                }
                glEnd();
                glDisable(GL_BLEND);
                glDisable(GL_POINT_SMOOTH);
            }
            
            // Calculate FPS
            frames++;
            uint32_t current_time = SDL_GetTicks();
            if (current_time - last_time >= 1000) {
                current_fps = frames * 1000.0f / (current_time - last_time);
                frames = 0;
                last_time = current_time;
            }

            // Stats Overlay HUD
            if (show_stats && ui_font) {
                char stats_buf[256];
                int peak_hz = (44100 / 2) / display_bands * peak_band;
                snprintf(stats_buf, sizeof(stats_buf), "FPS: %.1f | RMS: %.2f | PeakFreq: ~%dHz", current_fps, rms, peak_hz);
                SDL_Color text_col = {0, 255, 100, 255};
                render_text_ortho(stats_buf, 10, 10, ui_font, text_col, window_width, window_height);
                
                snprintf(stats_buf, sizeof(stats_buf), "Bass Eng: %.2f | [%s]", bass_sum, mode_names[current_mode_index] + 20);
                render_text_ortho(stats_buf, 10, 40, ui_font, text_col, window_width, window_height);
            }
            
            // Swap the double buffer
            SDL_GL_SwapWindow(window);
        }
        
        // Short delay ~60fps
        SDL_Delay(16);
    }

    // 6. Cleanup
    printf("\nPlayback and visualization finished.\n");
    cleanup_audio_player(&player);
    cleanup_gpu_dsp();
    free_media(&audio_data);
    free(magnitudes);
    free(chromagram);
    free(particles);
    free(galaxy_particles);
    free(grid_particles); // Added missing free for Mode 1
    
    for (int i=0; i<mountain_depth; i++) free(mountain_history[i]);
    free(mountain_history);

    if (ui_font) TTF_CloseFont(ui_font);
    TTF_Quit();
    gluDeleteQuadric(planet_quad);

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
