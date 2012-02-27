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
    //char host[] = "mdsplus.aug.ipp.mpg.de:8000";
    int sock;

    if ((sock=sm_mdsconnect(host)) < 0) {
	fprintf(stderr, "sock = %d\n", sock);
        sm_error("Could not connect.\n");
    }

    printf("sock = %d\n", sock);

    //sm_mdsopen(sock, "rcp", 132777);

    char cmd[] = "$";

    Descrip l, *R; 
	    
    R = (Descrip*) malloc(2*sizeof(Descrip));

    mkDescrip(&R[0], w_dtype_CSTRING, 0, NULL, 0, sizeof(char), cmd);

    int one = 1;
    double x = 1.;
    mkDescrip(&R[1], w_dtype_DOUBLE, 1, &one, 1, sizeof(double), &x);


    void *mem;
    sm_mdsvalue(sock, &l, 2, R, &mem);

    printf("out = %f\n", *(double*)l.ptr);

    if (mem) free(mem);
    free(R);

    //sm_mdsclose(sock);

    sm_mdsdisconnect(sock);

    exit(0);
}
