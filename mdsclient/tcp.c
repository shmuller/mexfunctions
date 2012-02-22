#ifdef _WIN32
#include "tcp_win32.h"
#else
#include "tcp_linux.h"
#endif

#include <stdlib.h>
#include <string.h>

#include <ipdesc.h>

int tcpopen(char *host, char *port) {   
    int sock;
    int one = 1;
    struct sockaddr_in sin;
    struct hostent *hp;

    if (Startup()) {
        return(-1);
    }
    if ((hp=gethostbyname(host)) == NULL) {
        return(-2);
    }

    memset(&sin,0,sizeof(sin));
    sin.sin_family=AF_INET;
    memcpy(&sin.sin_addr,hp->h_addr,hp->h_length);
    sin.sin_port = htons(atoi(port));

    if ((sock=socket(AF_INET,SOCK_STREAM,0)) < 0) {
        return(-3);
    }
    if (connect(sock,(struct sockaddr*)&sin,sizeof(sin)) < 0) {
        return(-4);
    }
    setsockopt(sock,IPPROTO_TCP,TCP_NODELAY,(void*)&one,sizeof(one));
    setsockopt(sock,SOL_SOCKET,SO_KEEPALIVE,(void*)&one,sizeof(one));
    setsockopt(sock,SOL_SOCKET,SO_OOBINLINE,(void*)&one,sizeof(one));
    return(sock);
}

int tcpauth(int sock) {
    int STRLEN = 4096;
    char user[STRLEN];
    char *user_p = getUserName(user,STRLEN);
    
    struct descrip exparg, *arg;
    int numbytes = 0, stat = 0;
    void *mem = NULL;
	
    arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,NULL,user_p);
    stat = SendArg(sock, 0, arg->dtype, 1, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
    if (!status_ok(stat)) {
        shutdown(sock,2);
        Cleanup();
        return(-5);
    }
    return(0);
}

int tcpconnect(char *host) {
    char *port;
    int sock, err;

    if ((port=strchr(host,':')) == NULL) {
        port = strdup("8000");
    } else {
        *port++ = 0;
    }
    if ((sock=tcpopen(host,port)) < 0) {
        return(sock);		
    }
    if ((err=tcpauth(sock)) < 0) {
        return(err);
    }
    return(sock);
}

int tcpdisconnect(int sock) {
    shutdown(sock,2);
    Cleanup();
    return(0);
}



