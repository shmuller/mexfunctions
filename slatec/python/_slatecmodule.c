#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>   /* memset */

extern void 
dintrv_(double *t, int *np1, double *x, int *inbv, int *i, int *mflag);


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
