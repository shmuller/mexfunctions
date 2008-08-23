#include <stdio.h>
#include "mdslib.h"
#include "ipdesc.h"

#include <winsock2.h>

int tcpopen(char *host, char *port)
{   
	int unit;
	struct sockaddr_in sin;
	struct hostent *hp;

	WSADATA wsadata;
	if (WSAStartup(MAKEWORD(1,1), &wsadata) == SOCKET_ERROR) {
		return(-1);
	}		
	
    	if ((hp=gethostbyname(host)) == NULL) {
		return(-2);
	}

	memset(&sin,0,sizeof(sin));
	sin.sin_family=AF_INET;
	memcpy(&sin.sin_addr,hp->h_addr,hp->h_length);
	sin.sin_port = htons(atoi(port));
	
	if ((unit=socket(AF_INET,SOCK_STREAM,0)) < 0) {
		return(-3);
	}
	if (connect(unit,(struct sockaddr*)&sin,sizeof(sin)) < 0) {
    		return(-4);
	}

	return(unit);
}


int main(int argc, char *argv[]) {

	int sock=0;
	if (argc > 2) {
		sock = tcpopen(argv[1],argv[2]);
	}
	printf("sock = %d\n",sock);

	if (sock > 0) {
		closesocket(sock);
	}

	WSACleanup();

/*
	SOCKET sock=0;

	if (argc > 1) {
		sock = ConnectToMds(argv[1]);
	}
	printf("sock = %d\n",(int)sock);

	if (sock == INVALID_SOCKET) {
		printf("Connection failed\n");
	}
*/

	return(0);
}
