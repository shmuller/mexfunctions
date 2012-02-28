#include <stdlib.h>
#include <string.h>

#include "tcp.h"

#ifdef _WIN32
#include "tcp_win32.h"
#else
#include "tcp_linux.h"
#endif

int tcpopen(char *host, char *port) {   
    int sock;
    int one = 1;
    struct sockaddr_in sin;
    struct hostent *hp;

    if (Startup()) {
        return -1;
    }
    if ((hp=gethostbyname(host)) == NULL) {
        return -2;
    }

    memset(&sin,0,sizeof(sin));
    sin.sin_family=AF_INET;
    memcpy(&sin.sin_addr,hp->h_addr,hp->h_length);
    sin.sin_port = htons(atoi(port));

    if ((sock=socket(AF_INET,SOCK_STREAM,0)) < 0) {
        return -3;
    }
    if (connect(sock,(struct sockaddr*)&sin,sizeof(sin)) < 0) {
        return -4;
    }
    setsockopt(sock,IPPROTO_TCP,TCP_NODELAY,(void*)&one,sizeof(one));
    setsockopt(sock,SOL_SOCKET,SO_KEEPALIVE,(void*)&one,sizeof(one));
    setsockopt(sock,SOL_SOCKET,SO_OOBINLINE,(void*)&one,sizeof(one));
    return sock;
}

int tcpclose(int sock) {
    shutdown(sock,2);
    Cleanup();
    return 0;
}



