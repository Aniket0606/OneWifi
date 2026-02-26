#ifndef MULTICAST_LISTENER_H
#define MULTICAST_LISTENER_H

#ifdef __cplusplus
extern "C" {
#endif

int multicast_listener_start(void);

void multicast_listener_stop(void);

int multicast_listener_is_running(void);

#ifdef __cplusplus
}
#endif

#endif /* MULTICAST_LISTENER_H */
