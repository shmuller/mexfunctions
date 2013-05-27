#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include "mdsclient.h"

int sm_error(const char *errstr) 
{
    fprintf(stderr, "%s\n", errstr);
    exit(0);
}

int main()
{
    char host[] = "localhost:8001";
    int sock;

    if ((sock=sm_mdsconnect(host)) < 0) {
	fprintf(stderr, "sock = %d\n", sock);
        sm_error("Could not connect.\n");
    }

    printf("sock = %d\n", sock);

    char cmd[] = "1+2";

    Descrip l, *R; 
	    
    R = (Descrip*) malloc(2*sizeof(Descrip));

    mkDescrip(&R[0], w_dtype_CSTRING, 0, NULL, 0, sizeof(char), cmd);

    //int one = 1;
    //double x = 1.;
    //mkDescrip(&R[1], w_dtype_DOUBLE, 1, &one, 1, sizeof(double), &x);


    sm_mdsvalue(sock, &l, 1, R);

    printf("out = %d\n", *(int*)l.ptr);

    if (l.ptr) free(l.ptr);
    free(R);

    sm_mdsdisconnect(sock);

    exit(0);
}
