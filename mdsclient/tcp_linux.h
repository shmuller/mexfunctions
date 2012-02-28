#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

#include <pwd.h>

int Startup() {return 0;}
int Cleanup() {return 0;}

char *tcpuser(char *user, int len)
{
    struct passwd pwd, *pwd_p;
    
    return getpwuid_r(geteuid(), &pwd, user, len, &pwd_p) ?
        "Linux User" : pwd_p->pw_name;
}

