#include <stdio.h>
#include "mdslib.h"
#include "ipdesc.h"

int main(int argc, char *argv[]) {

	SOCKET sock=0;

	if (argc > 1) {
		sock = ConnectToMds(argv[1]);
	}
	printf("sock = %d\n",(int)sock);

	if (sock == INVALID_SOCKET) {
		printf("Connection failed\n");
	}

	return(0);
}
