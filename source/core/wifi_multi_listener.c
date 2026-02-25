#include "multicast_listener.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <errno.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/time.h>

#include <cjson/cJSON.h>

#define MULTICAST_GROUP "225.0.0.50"
#define PORT 12575
#define BUFFER_SIZE 16 * 1024

static volatile int running = 0;
static pthread_mutex_t state_mutex = PTHREAD_MUTEX_INITIALIZER;

static void handle_json(char *json_str)
{
    cJSON *root = cJSON_Parse(json_str);
    if (!root)
        return;

    cJSON *_evCode = cJSON_GetObjectItem(root, "_evCode");

    if (_evCode && _evCode->valueint == 303) {

        cJSON *event = cJSON_GetObjectItem(
            root, "DeviceServiceResourceUpdatedEvent");

        if (event) {
            cJSON *resource = cJSON_GetObjectItem(event, "resource");
            if (resource) {
                cJSON *dsres = cJSON_GetObjectItem(resource, "DSResource");
                if (dsres) {
                    cJSON *id = cJSON_GetObjectItem(dsres, "id");
                    cJSON *value = cJSON_GetObjectItem(dsres, "value");

                    if (id && value &&
                        strcmp(id->valuestring, "faulted") == 0) {

                        wifi_util_info_print(WIFI_CTRL,"Door Sensor State: %s\n",
                               strcmp(value->valuestring, "true") == 0 ?
                               "OPEN (FAULT)" : "CLOSED (RESTORE)");
                    }
                }
            }
        }
    }

    cJSON_Delete(root);
}

static void* listener_thread_func(void *arg)
{
    int sockfd = -1;
    struct sockaddr_in addr;
    struct ip_mreq mreq;
    char buffer[BUFFER_SIZE];

    sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0) {
        perror("socket");
        goto exit_thread;
    }

    int reuse = 1;
    setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR,
               &reuse, sizeof(reuse));

    /* 1 second receive timeout */
    struct timeval tv;
    tv.tv_sec = 1;
    tv.tv_usec = 0;

    setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO,
               &tv, sizeof(tv));

    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(PORT);
    addr.sin_addr.s_addr = htonl(INADDR_ANY);

    if (bind(sockfd, (struct sockaddr*)&addr,
             sizeof(addr)) < 0) {
        perror("bind");
        goto exit_thread;
    }

    mreq.imr_multiaddr.s_addr = inet_addr(MULTICAST_GROUP);
    mreq.imr_interface.s_addr = htonl(INADDR_ANY);

    if (setsockopt(sockfd, IPPROTO_IP,
                   IP_ADD_MEMBERSHIP,
                   &mreq, sizeof(mreq)) < 0) {
        perror("IP_ADD_MEMBERSHIP");
        goto exit_thread;
    }

    wifi_util_info_print(WIFI_CTRL,"Detached multicast listener started\n");

    while (1) {

        pthread_mutex_lock(&state_mutex);
        int is_running = running;
        pthread_mutex_unlock(&state_mutex);

        if (!is_running)
            break;

        ssize_t len = recvfrom(sockfd,
                               buffer,
                               BUFFER_SIZE - 1,
                               0, NULL, NULL);

        if (len < 0) {
            if (errno == EAGAIN || errno == EWOULDBLOCK) {
                continue;  // timeout
            }
            perror("recvfrom");
            break;
        }

        buffer[len] = '\0';

        wifi_util_info_print(WIFI_CTRL,"====>recv message:%s\n", buffer);

        char *json_start = strchr(buffer, '|');
        if (!json_start)
            continue;

        json_start++;  // skip '|'
        while (*json_start == ' ')
            json_start++;

        handle_json(json_start);
    }

exit_thread:

    if (sockfd >= 0)
        close(sockfd);

    pthread_mutex_lock(&state_mutex);
    running = 0;
    pthread_mutex_unlock(&state_mutex);

    wifi_util_info_print(WIFI_CTRL,"Detached multicast listener exiting\n");

    pthread_exit(NULL);
}

int multicast_listener_start(void)
{
    pthread_mutex_lock(&state_mutex);

    if (running) {
        pthread_mutex_unlock(&state_mutex);
        return 0;  // already running
    }

    running = 1;

    pthread_t thread;
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

    if (pthread_create(&thread,
                       &attr,
                       listener_thread_func,
                       NULL) != 0) {
        running = 0;
        pthread_attr_destroy(&attr);
        pthread_mutex_unlock(&state_mutex);
        return -1;
    }

    pthread_attr_destroy(&attr);

    pthread_mutex_unlock(&state_mutex);
    return 0;
}

void multicast_listener_stop(void)
{
    pthread_mutex_lock(&state_mutex);
    running = 0;
    pthread_mutex_unlock(&state_mutex);
}

int multicast_listener_is_running(void)
{
    pthread_mutex_lock(&state_mutex);
    int state = running;
    pthread_mutex_unlock(&state_mutex);
    return state;
}
