# Set the audio backend. Options: ALSA, PULSEAUDIO
AUDIO_BACKEND ?= PULSEAUDIO

NVCC = nvcc
CFLAGS = -O3 -I/usr/include/SDL2
LDFLAGS_SDL = -lSDL2 -lavformat -lavcodec -lavutil -lswresample -lGL -lGLU -lSDL2_ttf

ifeq ($(AUDIO_BACKEND), PULSEAUDIO)
    CFLAGS += -DUSE_PULSEAUDIO
    LDFLAGS = $(LDFLAGS_SDL) -lpulse-simple -lpulse -lm
else
    CFLAGS += -DUSE_ALSA
    LDFLAGS = $(LDFLAGS_SDL) -lasound -lm
endif
SRC_DIR = src
OBJ_DIR = obj

# Find all C, CPP, and CU files in the source directory
C_SOURCES = $(wildcard $(SRC_DIR)/*.c)
CPP_SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
CU_SOURCES = $(wildcard $(SRC_DIR)/*.cu)

# Generate object file names
C_OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(C_SOURCES))
CPP_OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(CPP_SOURCES))
CU_OBJECTS = $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(CU_SOURCES))
OBJECTS = $(C_OBJECTS) $(CPP_OBJECTS) $(CU_OBJECTS)

TARGET = cuda_music_player

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
