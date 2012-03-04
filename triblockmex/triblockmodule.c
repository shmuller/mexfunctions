/* triblock(A,B,IPIV,C,D)
 * Tri-block-diagnonal solve step. Compile with:
 *
 * python setup.py install
 *
 * S. H. Muller, 2012/03/04
 */

#include <Python.h>
#include <numpy/arrayobject.h>


extern void dgesv(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
    double *B, int *LDB, int *INFO);

extern void dgetrf(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

extern void dgetrs(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
    int *IPIV, double *B, int *LDB, int *INFO);

extern void dgemm(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
    double *ALPHA, double *A, int *LDA, double *B, int *LDB,
    double *BETA, double *C, int *LDC);


static PyObject *triblock(PyObject *self, PyObject *args)
{
    PyObject *AA, *BB, *IP, *CC, *DD;
    if (!PyArg_ParseTuple(args, "OOOOO", &AA, &BB, &IP, &CC, &DD)) {
        PyErr_SetString(PyExc_TypeError, "A, B, IPIV, C, D expected");
        return NULL;
    }

	double *A = PyArray_DATA(AA);
	double *B = PyArray_DATA(BB);
	int *IPIV = PyArray_DATA(IP);
	double *C = PyArray_DATA(CC);
	double *D = PyArray_DATA(DD);

	int N = PyArray_DIM(BB,0), NRHS = PyArray_DIM(BB,1);

	int LDA = N, LDB = N, INFO = 0, N2 = N/2, NN2 = N*N2;
	char TRANSN = 'N', TRANST = 'T';
	double ALPHA = -1.0, BETA = 1.0;

	dgetrf(&N, &N, A, &LDA, IPIV, &INFO);
	
	if (INFO == 0) {
 		dgetrs(&TRANST, &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
 	}
        
 	D += NN2;
 	dgemm(&TRANST, &TRANSN, &N, &N2, &N, &ALPHA, B, &N, C, &N, &BETA, D, &N);
	
	Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
    {"triblock", triblock, METH_VARARGS, "Tri-block-diagonal solve step"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
inittriblock(void)
{
    Py_InitModule("triblock", methods);
}

