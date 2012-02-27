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
    char host[] = "localhost:8010";
    int sock;

    if ((sock=sm_mdsconnect(host)) < 0) {
        sm_error("Could not connect.\n");
    }

    printf("sock = %d\n", sock);

    sm_mdsopen(sock, "rcp", 132777);

    char cmd[] = "$shot";

    Descrip l, r, *R = &r;

    mkDescrip(R, w_dtype_CSTRING, 0, NULL, 0, sizeof(char), cmd);

    void *mem;
    sm_mdsvalue(sock, &l, 1, R, &mem);

    if (l.w_dtype == w_dtype_LONG)
        printf("out = %d\n", *(int*)l.ptr);

    if (mem) free(mem);
 
    sm_mdsclose(sock);

    sm_mdsdisconnect(sock);

    exit(0);
}
