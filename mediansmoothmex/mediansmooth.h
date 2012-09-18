#ifdef __cplusplus
extern "C"
{
#endif

int median_filt_double(const double *x, const int N, const int w, int *ind);
int median_filt_float(const float *x, const int N, const int w, int *ind);

#include "mediansmooth_pqueue.h"

#include <stdio.h>

// DOUBLE
static int lt_double(void *a, void *b) {
	return *(double*)a < *(double*)b;
}

static int gt_double(void *a, void *b) {
	return *(double*)a > *(double*)b;
}

static void print_double(void *a) {
    printf("%f", *(double*)a);
}

static double to_double_double(void *a) {
    return *(double*)a;
}

static void from_double_double(void *a, double x) {
    *(double*)a = x;
}

void median_filt_pqueue_double(void *X, int N, int w, int bdry) {
    fun_t fun = {&lt_double, &gt_double, &print_double, 
        &to_double_double, &from_double_double};
    median_filt_pqueue(X, N, w, bdry, sizeof(double), &fun);
}

// FLOAT
static int lt_float(void *a, void *b) {
	return *(float*)a < *(float*)b;
}

static int gt_float(void *a, void *b) {
	return *(float*)a > *(float*)b;
}

static void print_float(void *a) {
    printf("%f", *(float*)a);
}

static double to_double_float(void *a) {
    return *(float*)a;
}

static void from_double_float(void *a, double x) {
    *(float*)a = x;
}

void median_filt_pqueue_float(void *X, int N, int w, int bdry) {
    fun_t fun = {&lt_float, &gt_float, &print_float, 
        &to_double_float, &from_double_float};
    median_filt_pqueue(X, N, w, bdry, sizeof(float), &fun);
}

// INT64
static int lt_int64(void *a, void *b) {
	return *(long long*)a < *(long long*)b;
}

static int gt_int64(void *a, void *b) {
	return *(long long*)a > *(long long*)b;
}

static void print_int64(void *a) {
    printf("%lld", *(long long*)a);
}

static double to_double_int64(void *a) {
    return *(long long*)a;
}

static void from_double_int64(void *a, double x) {
    *(long long*)a = x;
}

void median_filt_pqueue_int64(void *X, int N, int w, int bdry) {
    fun_t fun = {&lt_int64, &gt_int64, &print_int64, 
        &to_double_int64, &from_double_int64};
    median_filt_pqueue(X, N, w, bdry, sizeof(long long), &fun);
}

// INT32
static int lt_int32(void *a, void *b) {
	return *(int*)a < *(int*)b;
}

static int gt_int32(void *a, void *b) {
	return *(int*)a > *(int*)b;
}

static void print_int32(void *a) {
    printf("%d", *(int*)a);
}

static double to_double_int32(void *a) {
    return *(int*)a;
}

static void from_double_int32(void *a, double x) {
    *(int*)a = x;
}

void median_filt_pqueue_int32(void *X, int N, int w, int bdry) {
    fun_t fun = {&lt_int32, &gt_int32, &print_int32, 
        &to_double_int32, &from_double_int32};
    median_filt_pqueue(X, N, w, bdry, sizeof(int), &fun);
}


#ifdef __cplusplus
}
#endif
