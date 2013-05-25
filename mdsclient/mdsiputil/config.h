/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated automatically from configure.in by autoheader.  */
/* @(#)$RCSfile: config.h.in,v $ $Revision: 1.15 $ */
#ifndef _CONFIG_INCLUDED

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#define TIME_WITH_SYS_TIME 1

/* Define if your <sys/time.h> declares struct tm.  */
/* #undef TM_IN_SYS_TIME */

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* _big_endian --- Compatibility with old definitions */
#ifdef WORDS_BIGENDIAN
#define _big_endian
#endif

#define TARGET_CPU i386
#define TARGET_VENDOR apple
#define TARGET_OS darwin11.2.0

/* Define NEED_SEMUN if sys/sem.h doesn't define union semun  */
/* #undef NEED_SEMUN */

/* Define NEED_SIGVEC if signal.h doesn't define struct sigvec  */
/* #undef NEED_SIGVEC */

/* Define PS_SYSV if ps(1) uses sysv arguments  */
#define PS_SYSV 1

/* Define TYPE_SHMAT if sys/shm.h doesn't type shmat()  */
/* #undef TYPE_SHMAT */

/* Define USE_NIS if system has nis_local_host() function  */
/* #undef USE_NIS */

/* Define PROTECTED if shared memory should be mapped read only*/
/* #undef PROTECTED */

/* Define NOCELL if IPCS is not part of an IPCS Cell*/
/* #undef NOCELL */

/* Define FILE_PTR_HL if the FILE structure has a high part and low part fd instead of just one  */
/* #undef FILE_PTR_HL */

/* The number of bytes in a long.  */
#define SIZEOF_LONG 8

/* The number of bytes in a pointer.  */
#define SIZEOF_INT_P 8

/* The number of bytes in a long long.  */
#define SIZEOF_LONG_LONG 8

/* The number of bytes in a _int64.  */
#define SIZEOF__INT64 0

/* Define if you have the gethostname function.  */
#define HAVE_GETHOSTNAME 1

/* Define if you have the gettimeofday function.  */
#define HAVE_GETTIMEOFDAY 1

/* Define if you have the select function.  */
#define HAVE_SELECT 1

/* Define if you have the sigvec function.  */
#define HAVE_SIGVEC 1

/* Define if you have the sigvector function.  */
/* #undef HAVE_SIGVECTOR */

/* Define if you have the socket function.  */
#define HAVE_SOCKET 1

/* Define if you have the <fcntl.h> header file.  */
#define HAVE_FCNTL_H 1

/* Define if you have the <netdb.h> header file.  */
#define HAVE_NETDB_H 1

/* Define if you have the <resolv.h> header file.  */
#define HAVE_RESOLV_H 1

/* Define if you have the <stdarg.h> header file.  */
#define HAVE_STDARG_H 1

/* Define if you have the <strings.h> header file.  */
#define HAVE_STRINGS_H 1

/* Define if you have the <sys/filio.h> header file.  */
#define HAVE_SYS_FILIO_H 1

/* Define if you have the <sys/ioctl.h> header file.  */
#define HAVE_SYS_IOCTL_H 1

/* Define if you have the <sys/resource.h> header file.  */
#define HAVE_SYS_RESOURCE_H 1

/* Define if you have sys <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define if you have the <syslog.h> header file.  */
#define HAVE_SYSLOG_H 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define if you have the <dl.h> header file. */
/* #undef HAVE_DL_H */

/* Define if you have the <vxWorks.h> header file. */
/* #undef HAVE_VXWORKS_H */

#if defined(HAVE_VXWORKS_H) && !defined(vxWorks)
#define vxWorks
#endif 

/* Set SHARELIB_TYPE to ".so",".sl",".a" etc... */
#define SHARELIB_TYPE ".dylib"

/* Define if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define if you have pthread_lock_global_np routine. */
/* #undef HAVE_PTHREAD_LOCK_GLOBAL_NP */

/* Define if you have pthread_mutexattr_settype routine. */
#define HAVE_PTHREAD_MUTEXATTR_SETTYPE 1

/* Define if you have pthread_mutexattr_setkind_np routine. */
/* #undef HAVE_PTHREAD_MUTEXATTR_SETKIND_NP */

/* Define if you have <readline/readline.h> header file. */
#define HAVE_READLINE_READLINE_H 1

/* Define if you have the getaddrinfo routine */
#define HAVE_GETADDRINFO 1

/* Is malloc.h available?  */
/* #undef HAVE_MALLOC_H */

/* Define if you prefer to use POSIX semaphores over SysV sem */
/* #undef USE_SEMAPHORE_H */

/* Define if you prefer to use named pipes over SysV msg */
/* #undef USE_PIPED_MESSAGING */

/* Define if you want to use new mdsip connections features */
/* #undef MDSIP_CONNECTIONS */

/* Define if you are using SRB */
/* #undef SRB */

/* Define if use need to use localtime instead of timezone and daylight */
#define USE_TM_GMTOFF 1

/* Define if you want to collect performance statistics -- linux only */
/* #undef USE_PERF */

/* Define if you compiler likes 0xnnnnll long long constants  */
#define HAVE_LL_CONSTANTS

#ifdef HAVE_LL_CONSTANTS
#define LONG_LONG_CONSTANT(value) value##ll
#else
#define LONG_LONG_CONSTANT(value) value
#endif

/* #undef HAS_SIGVEC */
#ifdef HAVE_SIGVEC
#define SIGVEC sigvec
#define HAS_SIGVEC
#endif
#ifdef HAVE_SIGVECTOR
#define SIGVEC sigvector
#define HAS_SIGVEC
#endif
#ifndef HAS_SIGVEC
#define SIGVEC sigvec
#endif
#define BLOCK_ALL 0xfffff

#if SIZEOF_INT_P == 8
#define _pointer_int _int64
#else
#define _pointer_int int
#endif

#ifdef HAVE_WINDOWS_H
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define _CONFIG_INCLUDED
#endif
