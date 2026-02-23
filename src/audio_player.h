#ifndef AUDIO_PLAYER_H
#define AUDIO_PLAYER_H

#include "media_reader.h"
#include <stdbool.h>
#include <stdint.h>
#include <pthread.h>

typedef struct {
    AudioData* audio_data;
    bool is_playing;
    uint32_t current_sample_index; // Track where we are in the playback
    void* backend_handle;
} AudioPlayer;

bool init_audio_player(AudioPlayer* player, AudioData* audio_data);
void start_playback(AudioPlayer* player);
void cleanup_audio_player(AudioPlayer* player);

#endif // AUDIO_PLAYER_H
