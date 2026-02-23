#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_ttf.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

extern "C" {
    #include "media_reader.h"
    #include "audio_player.h"
    #include "visualizer_types.h"
}
#include "gpu_dsp.h"

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
        fprintf(stderr, "OpenGL context could not be created! SDL Error: %s\n", SDL_GetError());
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
    
    SDL_GL_SetSwapInterval(1);
    
    glClearColor(20.0f/255.0f, 20.0f/255.0f, 25.0f/255.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // 5. Start Playback
    printf("Starting playback...\n");
    start_playback(&player);

    uint32_t num_frames_per_chunk = 8192;
    uint16_t num_channels = audio_data.header.num_channels;
    
    float* magnitudes = (float*)malloc(MAX_DISPLAY_BANDS * 2 * sizeof(float)); 
    uint32_t display_bands_returned;
    
    float* chromagram = (float*)malloc(12 * 2 * sizeof(float));
    uint32_t chroma_pitches;
    
    int mountain_depth = 100;
    float** mountain_history = (float**)malloc(mountain_depth * sizeof(float*));
    for (int i=0; i<mountain_depth; ++i) {
        mountain_history[i] = (float*)calloc(MAX_DISPLAY_BANDS * 2, sizeof(float));
    }
    
    // Setup for Volumetric Particle Cloud
    int num_particles = 10000;
    Particle* particles = (Particle*)calloc(num_particles, sizeof(Particle));
    for (int i = 0; i < num_particles; ++i) {
        particles[i].x = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].y = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].z = ((rand() % 2000) - 1000) / 10.0f;
        particles[i].life = (rand() % 100) / 100.0f;
    }

    // Setup for Ultimate Procedural Galaxy
    int num_galaxy_particles = 20000;
    Particle* galaxy_particles = (Particle*)calloc(num_galaxy_particles, sizeof(Particle));
    for (int i = 0; i < num_galaxy_particles; ++i) {
        float angle = ((float)rand() / RAND_MAX) * 20.0f * M_PI;
        float radius = ((float)rand() / RAND_MAX) * 250.0f + (angle * 4.0f);
        galaxy_particles[i].x = cosf(angle) * radius;
        galaxy_particles[i].y = ((float)rand() / RAND_MAX) * 20.0f - 10.0f;
        galaxy_particles[i].z = sinf(angle) * radius;
        galaxy_particles[i].vx = angle;
        galaxy_particles[i].vy = radius;
        galaxy_particles[i].vz = ((float)rand() / RAND_MAX) * M_PI * 2.0f;
    }

    // Setup for Straight Line Grid (Mode 1)
    int num_lines = 600;
    int points_per_line = 1200;
    int total_grid_particles = num_lines * points_per_line;
    Particle* grid_particles = (Particle*)calloc(total_grid_particles, sizeof(Particle));
    
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
            grid_particles[idx].vx = 0.0f;
            grid_particles[idx].vy = 0.0f;
            grid_particles[idx].vz = 0.0f;
            grid_particles[idx].life = 0.0f;
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

    // Zoom/Pan state for Mode 0
    float zoom_level = 1.0f;   // 1.0 = see all bands
    float pan_offset = 0.0f;   // 0.0 = leftmost
    int mouse_x = 0, mouse_y = 0;
    bool mouse_dragging = false;
    int drag_start_x = 0;
    float drag_start_pan = 0.0f;

    // Audio analysis globals (for HUD)
    float rms = 0;
    float bass_sum = 0;
    float peak_mag = 0;
    int peak_band = 0;

    while (player.is_playing && !quit) {
        if (current_mode_index != last_mode_index) {
            SDL_SetWindowTitle(window, mode_names[current_mode_index]);
            last_mode_index = current_mode_index;
        }

        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
                player.is_playing = false;
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
                else if (e.key.keysym.sym >= SDLK_1 && e.key.keysym.sym <= SDLK_9) {
                    current_mode_index = e.key.keysym.sym - SDLK_1;
                }
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
                // Playback Seeking: [ and ] keys
                else if (e.key.keysym.sym == SDLK_LEFTBRACKET) {
                    // Seek backward 5 seconds
                    uint32_t seek_frames = (uint32_t)(audio_data.header.sample_rate * 5);
                    if (player.current_sample_index > seek_frames) {
                        player.current_sample_index -= seek_frames;
                    } else {
                        player.current_sample_index = 0;
                    }
                }
                else if (e.key.keysym.sym == SDLK_RIGHTBRACKET) {
                    // Seek forward 5 seconds
                    uint32_t seek_frames = (uint32_t)(audio_data.header.sample_rate * 5);
                    player.current_sample_index += seek_frames;
                    if (player.current_sample_index >= audio_data.num_samples) {
                        player.current_sample_index = audio_data.num_samples - 1;
                    }
                }
            }
            // Mouse wheel zoom for Mode 0
            if (e.type == SDL_MOUSEWHEEL) {
                if (e.wheel.y > 0) {
                    zoom_level *= 1.2f;
                    if (zoom_level > 16.0f) zoom_level = 16.0f;
                } else if (e.wheel.y < 0) {
                    zoom_level /= 1.2f;
                    if (zoom_level < 1.0f) zoom_level = 1.0f;
                }
                // Clamp pan after zoom change
                float max_pan = 1.0f - (1.0f / zoom_level);
                if (pan_offset > max_pan) pan_offset = max_pan;
            }
            // Mouse drag for panning
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                mouse_dragging = true;
                drag_start_x = e.button.x;
                drag_start_pan = pan_offset;
            }
            if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                mouse_dragging = false;
            }
            if (e.type == SDL_MOUSEMOTION) {
                mouse_x = e.motion.x;
                mouse_y = e.motion.y;
                if (mouse_dragging && zoom_level > 1.0f) {
                    float delta = (float)(drag_start_x - e.motion.x) / (float)window_width;
                    pan_offset = drag_start_pan + delta;
                    float max_pan = 1.0f - (1.0f / zoom_level);
                    if (pan_offset < 0.0f) pan_offset = 0.0f;
                    if (pan_offset > max_pan) pan_offset = max_pan;
                }
            }
        }

        uint32_t play_index = player.current_sample_index;
        
        if (play_index + num_frames_per_chunk < audio_data.num_samples) {
            size_t sample_offset = play_index * num_channels;
            
            process_audio_chunk(&audio_data.audio_data[sample_offset], num_frames_per_chunk, num_channels, magnitudes, &display_bands_returned);
            int display_bands = display_bands_returned;
            
            // Shift mountain history
            for (int r = mountain_depth - 1; r > 0; --r) {
                memcpy(mountain_history[r], mountain_history[r-1], MAX_DISPLAY_BANDS * 2 * sizeof(float));
            }
            memcpy(mountain_history[0], magnitudes, MAX_DISPLAY_BANDS * 2 * sizeof(float));
            
            // Audio analysis for HUD
            rms = 0;
            bass_sum = 0;
            peak_mag = 0;
            peak_band = 0;
            for(int i = 0; i < display_bands; ++i) {
                float mag = (magnitudes[i * num_channels + 0] + magnitudes[i * num_channels + (num_channels>1?1:0)]) / 2.0f;
                rms += mag * mag;
                if(i < 5) bass_sum += mag;
                if(mag > peak_mag) { peak_mag = mag; peak_band = i; }
            }
            rms = sqrtf(rms / display_bands);

            // Clear the frame
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST);

            // Build the RenderContext for all modes
            RenderContext ctx;
            ctx.magnitudes = magnitudes;
            ctx.mountain_history = mountain_history;
            ctx.audio_raw = audio_data.audio_data;
            ctx.sample_offset = sample_offset;
            ctx.display_bands = display_bands;
            ctx.num_channels = num_channels;
            ctx.num_frames = num_frames_per_chunk;
            ctx.window_width = window_width;
            ctx.window_height = window_height;
            ctx.mountain_depth = mountain_depth;
            ctx.particles = particles;
            ctx.grid_particles = grid_particles;
            ctx.galaxy_particles = galaxy_particles;
            ctx.num_particles = num_particles;
            ctx.num_grid_particles = total_grid_particles;
            ctx.num_galaxy_particles = num_galaxy_particles;
            ctx.num_lines = num_lines;
            ctx.points_per_line = points_per_line;
            ctx.planet_quad = planet_quad;
            ctx.zoom_level = zoom_level;
            ctx.pan_offset = pan_offset;
            ctx.mouse_x = mouse_x;
            ctx.mouse_y = mouse_y;
            ctx.sample_rate = (float)audio_data.header.sample_rate;
            ctx.total_samples = audio_data.num_samples;
            ctx.current_sample = play_index;
            ctx.font = ui_font;

            // Dispatch to the correct mode's render function
            if (current_mode_index == MODE_BARS) {
                render_mode_bars(&ctx);
            } else if (current_mode_index == MODE_RADIAL) {
                render_mode_radial(&ctx);
            } else if (current_mode_index == MODE_CHROMAGRAM) {
                render_mode_stereo_pond(&ctx);
            } else if (current_mode_index == MODE_SPECTROGRAM) {
                render_mode_spectrogram(&ctx);
            } else if (current_mode_index == MODE_OSCILLOSCOPE) {
                render_mode_oscilloscope(&ctx);
            } else if (current_mode_index == MODE_PHASE_SPACE) {
                render_mode_phase_space(&ctx);
            } else if (current_mode_index == MODE_3D_MOUNTAIN) {
                render_mode_mountain(&ctx);
            } else if (current_mode_index == MODE_PARTICLE_CLOUD) {
                render_mode_particle_cloud(&ctx);
            } else if (current_mode_index == MODE_3D_LISSAJOUS) {
                render_mode_lissajous(&ctx);
            } else if (current_mode_index == MODE_ULTIMATE_GALAXY) {
                render_mode_galaxy(&ctx);
            }
            
            // FPS Calculation
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
            
            // Progress Bar (bottom of screen)
            {
                float progress = (float)play_index / (float)audio_data.num_samples;
                int bar_h = 4;
                int bar_y = window_height - bar_h;

                glMatrixMode(GL_PROJECTION); glLoadIdentity();
                gluOrtho2D(0, window_width, window_height, 0);
                glMatrixMode(GL_MODELVIEW); glLoadIdentity();

                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                // Background track
                glColor4ub(80, 80, 80, 150);
                glBegin(GL_QUADS);
                glVertex2f(0, bar_y);
                glVertex2f(window_width, bar_y);
                glVertex2f(window_width, bar_y + bar_h);
                glVertex2f(0, bar_y + bar_h);
                glEnd();

                // Progress fill
                glColor4ub(0, 220, 120, 220);
                glBegin(GL_QUADS);
                glVertex2f(0, bar_y);
                glVertex2f(progress * window_width, bar_y);
                glVertex2f(progress * window_width, bar_y + bar_h);
                glVertex2f(0, bar_y + bar_h);
                glEnd();

                glDisable(GL_BLEND);
            }
            
            SDL_GL_SwapWindow(window);
        }
        
        SDL_Delay(16);
    }

    // 6. Cleanup
    printf("\nPlayback and visualization finished.\n");
    player.is_playing = false;
    SDL_Delay(100); // Let PulseAudio thread drain before destroying resources
    cleanup_audio_player(&player);
    cleanup_gpu_dsp();
    free_media(&audio_data);
    free(magnitudes);
    free(chromagram);
    free(particles);
    free(galaxy_particles);
    free(grid_particles);
    
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
