/* mdsclient
 *
 * Compile on Linux: 
 * mex -v mdsclient.c COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 -I/usr/local/mdsplus/include -L/usr/local/mdsplus/lib -lMdsIpShr
 * 
 * Compile on Windows:
 * mex -v mdsclient.c OPTIMFLAGS=-O3 -I"C:\PROGRA~1\MDSplus\DEVTOOLS\include" -L"C:\PROGRA~1\MDSplus\DEVTOOLS\lib" -lMdsIpShr
 *
 * S. H. Muller, 2008/02/05
 */

#include <stdio.h>
#include <stdlib.h>

#include <mex.h>
#include <matrix.h>
#include <ipdesc.h>
#include <string.h>

#include <winsock2.h>

#ifndef min
#define min(a,b) ((a)<(b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a)>(b) ? (a) : (b))
#endif

#define status_ok(status) (((status) & 1) == 1)

#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif

#ifdef COMPAT_V70
typedef int mwSize;
#endif

void* getnumarg(const mxArray *r, mxClassID id) {
	if (mxGetClassID(r) == id) { 
		return mxGetData(r);
	} else {
		mexErrMsgTxt("Wrong argument type");
	}
}

void *getstringarg(const mxArray *r) {
	short len = mxGetNumberOfElements(r);
	void *str = malloc((len+1)*sizeof(char));

	mxGetString(r,str,len+1);
	return str;
}

void *mex2mds_cmplx(const mxArray *r) {
    size_t i,num,siz;
    void *pr,*pi,*buf,*b;
    
    num = mxGetNumberOfElements(r);
    siz = mxGetElementSize(r);
    pr  = mxGetData(r);
    pi  = mxGetImagData(r);
    buf = malloc(2*num*siz);
    
    for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
        memcpy(b    ,pr,siz);
        memcpy(b+siz,pi,siz);
    }
    return buf;
}

void mds2mex_cmplx(const mxArray *r, void *buf) {
    size_t i,num,siz;
    void *pr,*pi,*b;
    
    num = mxGetNumberOfElements(r);
    siz = mxGetElementSize(r);
    pr  = mxGetData(r);
    pi  = mxGetImagData(r);
    
    for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
        memcpy(pr,b    ,siz);
        memcpy(pi,b+siz,siz);
    }
}

void *getarg(const mxArray *r, char *ndims, int **dims, mxClassID *mxID, mxComplexity *mxCo) {
    int i;
    const mwSize *dimsR = mxGetDimensions(r);
    *ndims = mxGetNumberOfDimensions(r);
    *mxID  = mxGetClassID(r);
    *mxCo  = (mxIsComplex(r)) ? mxCOMPLEX : mxREAL;
    
    if (*mxID != mxCHAR_CLASS) {
        /* remove singleton dimensions */
        for(i=*ndims-1; i>=0; i--) if (dimsR[i]==1) (*ndims)--; else break;
        
        *dims  = calloc(*ndims,sizeof(int));
        for(i=0; i<*ndims; i++) (*dims)[i]=dimsR[i];
        
        if (*mxCo==mxREAL) {
            return mxGetData(r);
        } else {
            if (*mxID != mxDOUBLE_CLASS && *mxID != mxSINGLE_CLASS) {
                mexErrMsgTxt("Complex data must be single or double");
            }
            return mex2mds_cmplx(r);
        }
    } else {
        *ndims = 0; *dims = NULL;
        return getstringarg(r);
    }
}

int mex2mds_type(mxClassID mxID, mxComplexity mxCo, char *dtype)
{
  switch (mxID)
  {
	case mxCHAR_CLASS    :  *dtype = DTYPE_CSTRING;    break;
	case mxUINT8_CLASS   :  *dtype = DTYPE_UCHAR;      break;
	case mxINT8_CLASS    :  *dtype = DTYPE_CHAR;       break;
	case mxUINT16_CLASS  :  *dtype = DTYPE_USHORT;     break;
	case mxINT16_CLASS   :  *dtype = DTYPE_SHORT;      break;
	case mxUINT32_CLASS  :  *dtype = DTYPE_ULONG;      break;
	case mxINT32_CLASS   :  *dtype = DTYPE_LONG;       break;
	case mxUINT64_CLASS  :  *dtype = DTYPE_ULONGLONG;  break;
	case mxINT64_CLASS   :  *dtype = DTYPE_LONGLONG;   break;
	case mxSINGLE_CLASS  :  *dtype = (mxCo==mxREAL) ? DTYPE_FLOAT  : DTYPE_COMPLEX;        break;
	case mxDOUBLE_CLASS  :  *dtype = (mxCo==mxREAL) ? DTYPE_DOUBLE : DTYPE_COMPLEX_DOUBLE; break;
  }
  return(1);
}

int mds2mex_type(struct descrip *d, mxClassID *mxID, mxComplexity *mxCo)
{
  switch (d->dtype)
  {
	case DTYPE_CSTRING         :  *mxID = mxCHAR_CLASS;    *mxCo = mxREAL;     break;
	case DTYPE_UCHAR           :  *mxID = mxUINT8_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_CHAR            :  *mxID = mxINT8_CLASS;    *mxCo = mxREAL;     break;
	case DTYPE_USHORT          :  *mxID = mxUINT16_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_SHORT           :  *mxID = mxINT16_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_ULONG           :  *mxID = mxUINT32_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_LONG            :  *mxID = mxINT32_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_ULONGLONG       :  *mxID = mxUINT64_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_LONGLONG        :  *mxID = mxINT64_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_FLOAT           :  *mxID = mxSINGLE_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_DOUBLE          :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_COMPLEX         :  *mxID = mxSINGLE_CLASS;  *mxCo = mxCOMPLEX;  break;
	case DTYPE_COMPLEX_DOUBLE  :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxCOMPLEX;  break;
	default                    :  *mxID = mxUNKNOWN_CLASS; *mxCo = mxREAL;     break;
  }
  return(1);
}

int mds2mex_dims(struct descrip *d, char *ndims, mwSize **dims)
{
	int i;
	*ndims = max(d->ndims,2);
	*dims = malloc(*ndims*sizeof(mwSize));

	for(i=0; i<d->ndims; i++) (*dims)[i] = max(d->dims[i],1);
	for(; i<*ndims; i++) (*dims)[i] = 1;
	return(1);
}

int tcpopen(char *host, char *port) {   
	int sock;
    int one = 1;
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

int tcpauth(sock) {
    typedef struct _msghdr {
        int msglen;
        int status;
        short length;
        unsigned char nargs;
        unsigned char descriptor_idx;
        unsigned char message_id;
        unsigned char dtype;
        signed char client_type;
        unsigned char ndims;
        int dims[MAX_DIMS];
        int fill;
    } MsgHdr;
    
    typedef struct _mds_message {
        MsgHdr h;
        char bytes[1];
    } Message;
    
    int SENDCAPABILITIES = 0xf, CompressionLevel = 0, bytes_sent = 0;
    
    char user_p[] = "smueller";
    int numbytes;
    
    Message *m = calloc(1,sizeof(MsgHdr)+strlen(user_p));
    m->h.client_type = SENDCAPABILITIES;
    m->h.length = strlen(user_p);
    m->h.msglen = sizeof(MsgHdr)+m->h.length;
    m->h.dtype = DTYPE_CSTRING;
    m->h.status = CompressionLevel;
    m->h.ndims = 0;
    memcpy(m->bytes,user_p,m->h.length);
    
    bytes_sent = send(sock,(char*)m,m->h.msglen,0);
    printf("len = %d, bytes_sent = %d\n",m->h.msglen,bytes_sent);
    
}

int sm_mdsconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	char *serv, *port;
    SOCKET sock;
    
    char user_p[] = "mdsuser";
    struct descrip exparg, *arg;
    int numbytes = 0, stat = 0;
    void *mem = NULL;

/*  if ((sock=ConnectToMds(serv)) == INVALID_SOCKET) {
        mexErrMsgTxt("Could not connect to server");
    }
*/
    serv = getstringarg(R[1]);
	if ((port=strchr(serv,':')) == NULL) {
        port = strdup("8000");
    } else {
        *port++ = 0;
    }
    if ((sock=tcpopen(serv,port)) < 0) {
        mexErrMsgTxt("Could not connect to server");
    }
    
    //tcpauth(sock);
    arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,NULL,user_p);
    stat = SendArg(sock, 0, arg->dtype, 1, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
    
	L[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	*((int*)mxGetData(L[0])) = sock;
	return(1);
}

int sm_mdsopen(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	char *tree = getstringarg(R[2]);
	int shot = *((int*)getnumarg(R[3],mxINT32_CLASS));
	int stat = MdsOpen(sock,tree,shot);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not open tree");
	}
*/	return(stat);
}

int sm_mdsclose(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = MdsClose(sock);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not close tree");
	}
*/	return(stat);
}

int sm_mdsdisconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = DisconnectFromMds(sock);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not disconnect from server");
	}
*/	return(stat);
}

int sm_mdsvalue(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	void *mxArg;
    mxClassID mxID;
	mxComplexity mxCo;
	
	struct descrip exparg, *arg;
	int i = 0, numbytes = 0, stat = 0;
	void *mem = 0, *out = 0;

	char dtype, ndims;
	int *dims;
    mwSize *dimsL;

    for(i=2; i<nR; i++) {
        mxArg = getarg(R[i],&ndims,&dims,&mxID,&mxCo);
        stat = mex2mds_type(mxID,mxCo,&dtype);
        arg = MakeDescrip(&exparg,dtype,ndims,dims,mxArg);
        stat = SendArg(sock, i-2, arg->dtype, nR-2, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    }
	stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
	stat = mds2mex_type(arg,&mxID,&mxCo);

	if (mxID == mxUNKNOWN_CLASS) {
		L[0] = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
	} else if (mxID == mxCHAR_CLASS) {
		out = calloc(numbytes+1,sizeof(char));
		memcpy(out,arg->ptr,numbytes);
		L[0] = mxCreateString((char*)out);
	} else {
		stat = mds2mex_dims(arg,&ndims,&dimsL);
		L[0] = mxCreateNumericArray(ndims,dimsL,mxID,mxCo);
        if (mxCo==mxREAL) {
            memcpy(mxGetData(L[0]),arg->ptr,numbytes);
        } else {
            mds2mex_cmplx(L[0],arg->ptr);
        }   
	}
	if (mem) free(mem);
	return(1);
}



void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	char *cmd = getstringarg(R[0]);
	int stat;

	if (strcmp(cmd,"mdsvalue")==0) {
		stat = sm_mdsvalue(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsconnect")==0) {
		stat = sm_mdsconnect(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsopen")==0) {
		stat = sm_mdsopen(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsclose")==0) {
		stat = sm_mdsclose(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsdisconnect")==0) {
		stat = sm_mdsdisconnect(nL,L,nR,R);
	} else {
		mexErrMsgTxt("Unknown command");
	}

	if (!status_ok(stat)) {
		mexErrMsgTxt("Untrapped error occurred");
	}
}

