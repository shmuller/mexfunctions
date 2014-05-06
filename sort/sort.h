#ifdef __cplusplus
extern "C" {
#endif

/* in sort.c */
void qsort_wrap(int *a, int n);
void insert_sort(int *a, int n);
void select_sort(int *a, int n);
void numpy_quicksort(int *a, int n);
void numpy_mergesort(int *a, int n);
void merge_sort(int *a, int n);
void quick_sort(int *a, int n);
void quick_select(int *a, int n, int k);
unsigned int bitmap_sort(unsigned int *a, int n, unsigned int *res);

/* in sort.cpp */
void stdsort(int *a, int n);
void timsort(int *a, int n);

#ifdef __cplusplus
}
#endif
