#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>   /* memset */

extern void 
dintrv_(double *t, int *np1, double *x, int *inbv, int *i, int *mflag);

extern void
dbspvn_(double *t, int *jhigh, int *k, int *index, double *x, int *ileft, 
        double *vnikx, double *work, int *iwork);


void nd_dot_product(double *a, int *s, double *b, int *n, int nd, 
                    double f, double *res)
{
    int s0, n0, ndm1, i;
    double *bnext, tmp;
    n0 = *n++;
    if (nd > 1) {
        s0 = *s++;
        ndm1 = nd-1;
        bnext = b+n0;
        for (i=n0+1; --i; a+=s0) {
            nd_dot_product(a, s, bnext, n, ndm1, f*(*b++), res);
        }
    } else {
        for (i=n0+1,tmp=0.; --i; tmp+=(*a++)*(*b++));
        *res += f*tmp;
    }
}


static PyObject* dbualnd(PyObject *self, PyObject *args)
{
    int ndim, *n, *k, *s, ideriv, m, *inbv, i, d, mm,
        jt, jb, jw, offs, nd, kd, np1, iwork, mflag, one=1;
    double *t, *a, *x, *work, *y, f;
    PyObject *obj;
    obj = PyTuple_GET_ITEM(args, 0);
    t = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 1);
    a = PyArray_DATA(obj);
   
    obj = PyTuple_GET_ITEM(args, 2);
    n = PyArray_DATA(obj);
    ndim = PyArray_SIZE(obj);

    obj = PyTuple_GET_ITEM(args, 3);
    k = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 4);
    s = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 5);
    ideriv = PyInt_AS_LONG(obj);

    obj = PyTuple_GET_ITEM(args, 6);
    x = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 7);
    inbv = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 8);
    work = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 9);
    y = PyArray_DATA(obj);
    m = PyArray_SIZE(obj);
    
    for (d=ndim,jw=k[d-1],s[d-1]=1; --d; jw+=k[d-1],s[d-1]=n[d]*s[d]);

    for (mm=0; mm<m; ++mm) {
        jt = jb = offs = 0;
        for (d=0; d<ndim; ++d) {
            nd = n[d];
            kd = k[d];
            np1 = nd+1;
            if (mm && x[0]==x[-ndim]) {
                i = inbv[d];
            } else {
                dintrv_(t+jt, &np1, x, inbv+d, &i, &mflag);
                if (mflag) i = (mflag == 1) ? nd : kd;
                dbspvn_(t+jt, &kd, &kd, &one, x, &i, work+jb, work+jw, &iwork);
            }
            jt += nd + kd;
            jb += kd;
            offs += (i-kd)*s[d];
            ++x;
        }
        f = 1.;
        *y = 0.;
        nd_dot_product(a+offs, s, work, k, ndim, f, y);
        ++y;
    }
    Py_RETURN_NONE;
}


/*
      SUBROUTINE DBUALND (NDIM, T, A, N, K, S, IDERIV, X, M, INBV, &
       WORK, Y)
      INTEGER NDIM, N(NDIM), K(NDIM), S(NDIM), INBV(*), I, OFFS, &
       IDERIV, KMIDER, IWORK, D, ND, KD, JT, JB, JW, M, MM
      DOUBLE PRECISION T(*), A(*), WORK(*), X(NDIM,M), Y(M), F
!***FIRST EXECUTABLE STATEMENT  DBUAL
      KMIDER = K(1) - IDERIV
      IF (KMIDER.LE.0) GO TO 99

      JW = SUM(K) + 1

      S(NDIM) = 1
      DO 10 D=NDIM,2,-1
        S(D-1) = N(D)*S(D)
   10 CONTINUE

      DO 50 MM=1,M
        JT = 1
        JB = 1
        OFFS = 1
        DO 20 D=1,NDIM
          ND = N(D)
          KD = K(D)
          IF (MM.GT.1 .AND. X(D,MM-1).EQ.X(D,MM)) THEN
            I = INBV(D)
          ELSE
            CALL DFINDI (T(JT), ND, KD, X(D,MM), INBV(D), I)
            CALL DBSPVN (T(JT), KD, KD, 1, X(D,MM), I, &
                         WORK(JB), WORK(JW), IWORK)
          END IF
          JT = JT + ND + KD
          JB = JB + KD
          OFFS = OFFS + (I-KD)*S(D)
   20   CONTINUE

        F = 1
        Y(MM) = 0
        CALL nd_dot_product(A(OFFS), S, WORK, K, NDIM, F, Y(MM))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 MM=1,M
  100   Y(MM) = 0.0D0
      RETURN
      END
*/

static PyObject* dbvali(PyObject *self, PyObject *args)
{
    int n, k, ideriv, m, *inbv, np1, i, mflag, mm, kmider, kk, j;
    double *t, *a, *x, *work, *y, *ti, *ai, f1, f2, fkmj, xm;
    PyObject *obj;
    obj = PyTuple_GET_ITEM(args, 0);
    t = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 1);
    a = PyArray_DATA(obj);
    n = PyArray_SIZE(obj);
   
    obj = PyTuple_GET_ITEM(args, 2);
    k = PyInt_AS_LONG(obj);

    obj = PyTuple_GET_ITEM(args, 3);
    ideriv = PyInt_AS_LONG(obj);

    obj = PyTuple_GET_ITEM(args, 4);
    x = PyArray_DATA(obj);
    m = PyArray_SIZE(obj);

    obj = PyTuple_GET_ITEM(args, 5);
    inbv = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 6);
    work = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 7);
    y = PyArray_DATA(obj);

    kmider = k - ideriv;
    if (kmider <= 0) {
        memset(y, 0., m);
        Py_RETURN_NONE;
    }
    np1 = n + 1;
    for (mm=m+1; --mm; ) {
        xm = *x++;
        dintrv_(t, &np1, &xm, inbv, &i, &mflag);
        if (mflag) i = (mflag == 1) ? n : k;
        ti = t + i;
        ai = a + i - k;

        for (kk=0; kk<k; ++kk) {
            work[kk] = ai[kk];
        }
        for (--kk; kk>=kmider; --kk) {
            fkmj = kk;
            for (j=0; j<kk; ++j) {
                work[j] = (work[j+1]-work[j])/(ti[j]-ti[j-kk])*fkmj;
            }
        }
        for (; kk>=1; --kk) {
            for(j=0; j<kk; ++j) {
                f1 = ti[j] - xm;
                f2 = ti[j-kk] - xm;
                work[j] = (work[j]*f1-work[j+1]*f2)/(f1-f2);
            }
        }
        *y++ = work[0];
    }
    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"dbualnd", dbualnd, METH_VARARGS, "Evaluate ND spline"},
    {"dbvali", dbvali, METH_VARARGS, "Evaluate B-spline"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_slatec(void)
{
    Py_InitModule("_slatec", methods);
    import_array();
}

/*
      SUBROUTINE DBVALI (T, A, N, K, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), K, KMIDER, KM1, N, NP1, &
       M, MM, MFLAG, KK, J, IMK, IPJ
      DOUBLE PRECISION T(*), A(N), WORK(*), X(*), Y(M), F1, F2, FKMJ
!***FIRST EXECUTABLE STATEMENT  DBVAL1
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KM1 = K - 1
      NP1 = N + 1
      DO 50 MM=1,M
        CALL DINTRV(T, NP1, X(MM), INBV(1), I, MFLAG)
        IF (MFLAG.NE.0) THEN
          IF (MFLAG.EQ.1) THEN
            I = N
          ELSE
            I = K
          END IF
        END IF
        IMK = I - K
        DO 5 KK=1,K
          WORK(KK) = A(IMK+KK)
    5   CONTINUE
        DO 20 KK=KM1,KMIDER,-1
          FKMJ = KK
          DO 10 J=1,KK
            IPJ = I + J
            WORK(J) = (WORK(J+1)-WORK(J))/(T(IPJ)-T(IPJ-KK))*FKMJ
   10     CONTINUE
   20   CONTINUE
        DO 40 KK=KMIDER-1,1,-1
          DO 30 J=1,KK
            IPJ = I + J
            F1 = T(IPJ) - X(MM)
            F2 = T(IPJ-KK) - X(MM)
            WORK(J) = (WORK(J)*F1-WORK(J+1)*F2)/(F1-F2)
   30     CONTINUE
   40   CONTINUE
        Y(MM) = WORK(1)
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,M
  100   Y(I) = 0.0D0
      RETURN
      END
*/

