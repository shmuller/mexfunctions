#ifndef __TCP_H__
#define __TCP_H__

#ifdef __cplusplus
extern "C" {
#endif

char *tcpuser(char *user, int len);

int tcpopen(char *host, char *port);
int tcpclose(int sock);

#ifdef __cplusplus
}
#endif

#endif /* __TCP_H__ */

