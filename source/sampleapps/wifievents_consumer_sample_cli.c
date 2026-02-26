#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "wifievents_consumer_sample.h"

// Map gesture ID to bitmask
uint32_t gesture_id_to_mask(int id)
{
    switch(id) {
        case 1: return GESTURE_HAND_MOVEMENT;
        case 2: return GESTURE_WALKING;
        case 3: return GESTURE_RUNNING;
        default: return 0;
    }
}

// Convert enum to string
const char* gesture_enum_to_string(motion_gesture_t gesture)
{
    switch (gesture) {
        case GESTURE_IDLE:         return "idle";
        case GESTURE_HAND_MOVEMENT: return "hand movement";
        case GESTURE_WALKING:       return "walking";
        case GESTURE_RUNNING:       return "running";
        // Add more gestures here
        default:                     return "unknown";
    }
}

void help(void)
{
    printf("Commands:\n");
    printf("  -g 0           -> stop all gestures\n");
    printf("  -g <id> <0|1> -> stop/start gesture\n");
    printf("Gesture IDs: 1=hand, 2=walking, 3=running\n");
}

// Print current gestures
void print_current_gestures(motion_gesture_t l_motion_gesture)
{
    printf("Current gestures: ");
    if (l_motion_gesture == 0) {
        printf("idle ");
    } else {
        if (l_motion_gesture & GESTURE_HAND_MOVEMENT) printf("hand ");
        if (l_motion_gesture & GESTURE_WALKING) printf("walking ");
        if (l_motion_gesture & GESTURE_RUNNING) printf("running ");
        if (l_motion_gesture & GESTURE_IDLE) printf("idle ");
    }
    printf("\n");
}

// Parse a single command line from serial
void process_serial_command(char *line)
{
    int g_id, state;
    motion_gesture_t l_motion_gesture = get_motion_gesture_obj();

    // Trim newline
    line[strcspn(line, "\n")] = 0;

    if (sscanf(line, "-g %d %d", &g_id, &state) == 2) {
        uint32_t mask = gesture_id_to_mask(g_id);
        if (mask == 0) {
            printf("Unknown gesture ID %d\n", g_id);
            help();
            return;
        }

        if (state == 1) {
            l_motion_gesture |= mask; // start gesture
            set_motion_gesture_obj(l_motion_gesture);
            printf("Gesture ID:%d[%s] started\n", g_id, gesture_enum_to_string(mask));
        } else if (state == 0) {
            l_motion_gesture &= ~mask; // stop gesture
            set_motion_gesture_obj(l_motion_gesture);
            printf("Gesture ID:%d[%s] stopped\n", g_id, gesture_enum_to_string(mask));
        } else {
            printf("Invalid state: %d\n", state);
            help();
        }
    } else if (sscanf(line, "-g %d", &g_id) == 1 && g_id == 0) {
        l_motion_gesture = 0; // stop all
        set_motion_gesture_obj(l_motion_gesture);
        printf("All gestures stopped\n");
    } else {
        printf("Invalid command: %s\n", line);
        printf("Use: -g 0 | -g <id> <0|1>\n");
        help();
    }

    print_current_gestures(l_motion_gesture);
}

void* serial_thread(void* arg) {
    char input_line[64];

    printf("Serial gesture input thread started.\n");
    help();

    while (1) {
        printf("> ");
        if (fgets(input_line, sizeof(input_line), stdin) != NULL) {
            process_serial_command(input_line);
        }
    }

    return NULL;
}

int start_cli_thread(void)
{
    pthread_t tid;

    // Start serial input thread
    if (pthread_create(&tid, NULL, serial_thread, NULL) != 0) {
        perror("Failed to create serial thread");
        return -1;
    } else {
        printf("cli thread started\r\n");
        pthread_detach(tid);
    }

    return 0;
}

