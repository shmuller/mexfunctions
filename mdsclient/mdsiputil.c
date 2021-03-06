#include "mdsip.h"
#ifdef BUFSIZ
#undef BUFSIZ
#define BUFSIZ 65536
#endif
#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#ifndef HAVE_VXWORKS_H
#include <pwd.h>
#endif
#endif
static unsigned char message_id = 1;
#ifdef NOCOMPRESSION
#define compress2(a,b,c,d,e) -1
#define uncompress(a,b,c,d) -1
static int CompressionLevel = 0;
#else
#include <zlib.h>
static int CompressionLevel = 0;
#endif

#define SocketSend send
#define SocketRecv recv
#define CloseSocket(s) -1
#define FlushSocket(s)

int SendBytes(SOCKET sock, char *bptr, int bytes_to_send, int options)
{
  int tries = 0;
  while ((bytes_to_send > 0) && (tries < 10))
  {
	int bytes_sent;
    int bytes_this_time = min(bytes_to_send,BUFSIZ);
    bytes_sent = SocketSend(sock, bptr, bytes_to_send, options);
    if (bytes_sent <= 0)
    {
      if (errno != EINTR)
        return 0;
      tries++;
    }
    else
    {
      bytes_to_send -= bytes_sent;
      bptr += bytes_sent;
      tries = 0;
    }
  }
  if (tries >= 10)
  {
    CloseSocket(sock);
    fprintf(stderr,"\rSendBytes shutdown socket %d",sock);
    return 0;
  }
  return 1;
}

int GetBytes(SOCKET sock, char *bptr, int bytes_to_recv, int oob)
{
  int tries = 0;
  while (bytes_to_recv > 0 && (tries < 10))
  {
    int bytes_recv;
    int bytes_this_time = min(bytes_to_recv,BUFSIZ);
    bytes_recv = SocketRecv(sock, bptr, bytes_to_recv, oob ? MSG_OOB : 0);
    if (bytes_recv <= 0)
    {
      if (errno != EINTR)
        return 0;
      tries++;
    }
    else
    {
      tries = 0;
      bytes_to_recv -= bytes_recv;
      bptr += bytes_recv;
    }
  }
  if (tries >= 10)
  {
    CloseSocket(sock);
    fprintf(stderr,"\rGetBytes shutdown socket %d: too many EINTR's",sock);
    return 0;
  }
  return 1;
}


int SetCompressionLevel(int level)
{
  int old_level = CompressionLevel;
  CompressionLevel = level;
  return old_level;
}

int GetCompressionLevel()
{
  return CompressionLevel;
}

char ClientType(void)
{
  static char ctype = 0;
  if (!ctype)
  {
    union { int i bits32;
            char  c[sizeof(double)];
            float x;
            double d;
          } client_test;
    client_test.x = 1.;
    if (client_test.i == 0x4080) {
      client_test.d=12345678;
      if(client_test.c[5])
       ctype = VMSG_CLIENT;
      else
      ctype = VMS_CLIENT;
    }
    else if (client_test.i == 0x3F800000)
    {
      if (sizeof(int) == 8)
        ctype = CRAY_IEEE_CLIENT;
      else
        ctype = IEEE_CLIENT;
    }
    else
      ctype = CRAY_CLIENT;
    client_test.i = 1;
    if (!client_test.c[0]) ctype |= BigEndian;
  }
  return ctype;
}

struct descrip *MakeDescrip(struct descrip *in_descrip, char dtype, char ndims, int *dims, void *ptr)
{
  int i;
  in_descrip->dtype = dtype;
  in_descrip->ndims = ndims;
  in_descrip->length = 0;
  for (i=0;i<ndims;i++) in_descrip->dims[i] = dims[i];
  for (i=ndims; i<MAX_DIMS; i++) in_descrip->dims[i] = 0;
  in_descrip->ptr = ptr;
  return in_descrip;
}

short ArgLen(struct descrip *d)
{
  short len;
  switch (d->dtype)
  {
    case DTYPE_CSTRING :  len = d->length ? d->length : (d->ptr ? strlen(d->ptr) : 0); break;
    case DTYPE_UCHAR   :
    case DTYPE_CHAR    :  len = sizeof(char); break;
    case DTYPE_USHORT  :
    case DTYPE_SHORT   :  len = sizeof(short); break;
    case DTYPE_ULONG   :  
    case DTYPE_LONG    :  len = sizeof(int); break;
    case DTYPE_FLOAT   :  len = sizeof(float); break;
    case DTYPE_DOUBLE  :  len = sizeof(double); break;
    case DTYPE_COMPLEX :  len = sizeof(float) * 2; break;
    case DTYPE_COMPLEX_DOUBLE :  len = sizeof(double) * 2; break;
    case DTYPE_ULONGLONG :
    case DTYPE_LONGLONG  :  len = 8; break;
  }
  return len;
}


#define FlipBytes(num,ptr) \
{\
  int __i;\
  int __n = num;\
  char *__p = ptr;\
  for (__i=0;__i<__n/2;__i++)\
  {\
    char __tmp = __p[__i];\
    __p[__i] = __p[__n - __i -1];\
    __p[__n - __i - 1] = __tmp;\
  }\
}

void FlipHeader(MsgHdr *header)
{
  int i;
#ifdef __CRAY
#define Flip32(n) n = ((n >> 24) & 0xff) | ((n >> 8) & 0xff00) | ((n << 8) & 0xff0000) | ((n << 24) & 0xff000000)
#define Flip16(n) n = ((n >> 8) & 0xff) | ((n << 8) & 0xff00)
  Flip32(header->msglen);
  Flip32(header->status);
  Flip16(header->length);
  for (i=0;i<MAX_DIMS;i++) FlipBytes(4,((char *)header->dims)+4*i);
#else
  FlipBytes(4,(char *)&header->msglen);
  FlipBytes(4,(char *)&header->status);
  FlipBytes(2,(char *)&header->length);
  for (i=0;i<MAX_DIMS;i++) FlipBytes(sizeof(header->dims[i]),(char *)&header->dims[i]);
#endif
}

static void FlipData(MsgHdr *h, char *bytes)
{
  int num = 1;
  int i;
  char *ptr;
  int dims[MAX_DIMS];
  for (i=0;i<MAX_DIMS;i++)
  {
#ifdef __CRAY
    dims[i] = i % 2 ? h->dims[i/2] & 0xffffffff : h->dims[i/2] >> 32;
#else
    dims[i] = h->dims[i];
#endif
  }
  if (h->ndims) for (i=0;i<h->ndims;i++) num *= dims[i];
#ifdef DEBUG
  printf("num to flip = %d\n",num);
#endif
  switch (h->dtype)
  {
#ifndef __CRAY
    case DTYPE_COMPLEX:
    case DTYPE_COMPLEX_DOUBLE: for (i=0,ptr=bytes;i<(num * 2);i++,ptr += h->length/2) FlipBytes(h->length/2,ptr); break;
    case DTYPE_FLOAT:   
    case DTYPE_DOUBLE:
#endif
    case DTYPE_LONGLONG:
    case DTYPE_ULONGLONG:
    case DTYPE_USHORT:
    case DTYPE_SHORT:  
    case DTYPE_ULONG:
    case DTYPE_LONG:       for (i=0,ptr=bytes;i<num;i++,ptr += h->length) FlipBytes(h->length,ptr); break;
  }
}


int SendMdsMsg(SOCKET sock, Message *m, int oob)
{
  unsigned long len = m->h.msglen - sizeof(m->h);
  unsigned long clength = 0; 
  Message *cm = 0;
  int status;
  int do_swap = 0; /*Added to handle byte swapping with compression*/
 
  if (len > 0 && CompressionLevel > 0 && m->h.client_type != SENDCAPABILITIES)
  {
	  clength = len;
	  cm = (Message *)malloc(m->h.msglen + 4);
  }
  if (!oob) FlushSocket(sock);
  if (m->h.client_type == SENDCAPABILITIES)
    m->h.status = CompressionLevel;
  if ((m->h.client_type & SwapEndianOnServer) != 0)
  {
    if (Endian(m->h.client_type) != Endian(ClientType()))
    {
      FlipData(&m->h, m->bytes);
      FlipHeader(&m->h);
      do_swap = 1; /* Recall that the header field msglen needs to be swapped */
    }
  }
  else
    m->h.client_type = ClientType();
  if (clength && compress2(cm->bytes+4,&clength,m->bytes,len,CompressionLevel) == 0 && clength < len)
  {
    cm->h = m->h;
    cm->h.client_type |= COMPRESSED;
    memcpy(cm->bytes,&cm->h.msglen,4);
    cm->h.msglen = clength + 4 + sizeof(MsgHdr);
/*If byte swapping required, swap msglen */
    if(do_swap)
 	FlipBytes(4,(char *)&cm->h.msglen);	
 /* status = SendBytes(sock, (char *)cm, cm->h.msglen, oob);*/
/* now msglen is swapped, and cannot be used as byte counter */
    status = SendBytes(sock, (char *)cm, clength + 4 + sizeof(MsgHdr), oob);
  }
  else
    status = SendBytes(sock, (char *)m, len + sizeof(MsgHdr), oob);
  if (clength) free(cm);
  return status;
}

int GetMdsMsg(SOCKET sock, MsgHdr *h, char **bytes)
{
  int msglen = 0;
  int status = GetBytes(sock, (char *)h, sizeof(MsgHdr), 0);
  if (status &1)
  {
    if ( Endian(h->client_type) != Endian(ClientType()) ) FlipHeader(h);
#ifdef DEBUG
    printf("msglen = %d\nstatus = %d\nlength = %d\nnargs = %d\ndescriptor_idx = %d\nmessage_id = %d\ndtype = %d\n",
               h->msglen,h->status,h->length,h->nargs,h->descriptor_idx,h->message_id,h->dtype);
    printf("client_type = %d\nndims = %d\n",h->client_type,h->ndims);
#endif
    if (CType(h->client_type) > CRAY_CLIENT || h->ndims > MAX_DIMS)
    {
      CloseSocket(sock);
      fprintf(stderr,"\rGetMdsMsg shutdown socket %d: bad msg header, header.ndims=%d, client_type=%d\n",sock,h->ndims,CType(h->client_type));
      *bytes = 0;
      return 0;
    }  
    msglen = h->msglen;
    *bytes = malloc(msglen - sizeof(MsgHdr));
    status = GetBytes(sock, *bytes, msglen - sizeof(MsgHdr), 0);
    
    if (status & 1 && IsCompressed(h->client_type))
    {
      char *b;
      unsigned long dlen;
      memcpy(&msglen, *bytes, 4);
      if (Endian(h->client_type) != Endian(ClientType()))
        FlipBytes(4,(char *)&msglen);
      dlen = msglen - sizeof(MsgHdr);
      b = malloc(dlen);
      status = uncompress(b, &dlen, *bytes + 4, h->msglen - sizeof(MsgHdr) - 4) == 0;
      if (status & 1)
      {
        h->msglen = msglen;
        free(*bytes);
        *bytes = b;
      }
      else
	    free(b);
    }
    
    if (status & 1 && (Endian(h->client_type) != Endian(ClientType())))
      FlipData(h, *bytes);
  }
  return status;
}


int  SendArg(SOCKET sock, unsigned char idx, char dtype, unsigned char nargs, short length, char ndims,
int *dims, char *bytes)
{
  int status;
  int msglen;
  int i;
  int nbytes = length;
  Message *m;
  if (idx > nargs)
  {
    /**** Special I/O message ****/ 
    nbytes = dims[0];
  }
  else
  {
    for (i=0;i<ndims;i++)
      nbytes *= dims[i];
  }
  msglen = sizeof(MsgHdr) + nbytes;
  m = memset(malloc(msglen),0,msglen);
  m->h.client_type = 0;
  m->h.msglen = msglen;
  m->h.message_id = message_id;
  m->h.descriptor_idx = idx;
  m->h.dtype = dtype;
  m->h.nargs = nargs;
  m->h.length = length;
  m->h.ndims = ndims;
#ifdef __CRAY
  for (i=0;i<4;i++)
    m->h.dims[i] = ((ndims > i * 2) ? (dims[i * 2] << 32) : 0) | ((ndims > (i * 2 + 1)) ? (dims[i * 2 + 1]) : 0); 
#else
  for (i=0;i<MAX_DIMS;i++) m->h.dims[i] = i < ndims ? dims[i] : 0;
#endif
  memcpy(m->bytes,bytes,nbytes);
  status = SendMdsMsg(sock, m, 0);
  if (idx > nargs || idx == (nargs -1)) message_id++;
  if (message_id == 0) message_id = 1;
  free(m);
  return status;
}

int  GetAnswerData(SOCKET sock, char *dtype, short *length, char *ndims, int *dims, int *numbytes, void **dptr)
{
  MsgHdr h;
  int i;
  *numbytes = 0;
  *dptr = 0;
  int status = GetMdsMsg(sock, &h, (char**) dptr);
  if (status != 1)
  {
    *dtype = 0;
    *length = 0;
    *ndims = 0;
    *numbytes = 0;
    if (*dptr)
    {
      free(*dptr);
      *dptr = 0;
    }
    return 0;
  }
  if (h.ndims)
  {
    *numbytes = h.length;
    for (i=0;i<h.ndims;i++)
    {
#ifdef __CRAY
      dims[i] = i % 2 ? h.dims[i/2] & 0xffffffff : h.dims[i/2] >> 32;
#else
      dims[i] = h.dims[i];
#endif
      *numbytes *= dims[i];
#ifdef DEBUG
      printf("dim[%d] = %d\n",i,dims[i]);
#endif
    }
    for (i=h.ndims;i < MAX_DIMS; i++)
      dims[i] = 0;
  }
  else
  {
    *numbytes = h.length;
    for (i=0;i<MAX_DIMS;i++)
      dims[i] = 0;
  }
  if ((int)(sizeof(MsgHdr) + *numbytes) != h.msglen)
  {
    *numbytes = 0;
    if (*dptr) {
      free(*dptr);
      *dptr=0;
    }
    return 0;
  }
  *dtype = h.dtype;
  *length = h.length;
  *ndims = h.ndims;
  return h.status;
}


