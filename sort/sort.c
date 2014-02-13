#include <stdio.h>
#include <stdlib.h> /* qsort */

#include "sort.h"

#define PYA_QS_STACK 100
#define SMALL_QUICKSORT 15
#define SMALL_MERGESORT 20

#define DTYPE_SWAP(a, b) {dtype tmp=(b); (b)=(a); (a)=tmp;}
#define DTYPE_LT(a, b) ((a) < (b))

typedef int dtype;

int compare_dtype(const void * a, const void * b) {
    return *(dtype*)a - *(dtype*)b;
}

void qsort_wrap(dtype *a, int n) {
    qsort(a, n, sizeof(dtype), compare_dtype);
}


void insert_sort(dtype *pl, int n) {
    dtype *pr = pl + n, *pi, *pj, *pk, vp;
    for (pi=pl+1; pi<pr; ++pi) {
        vp = *pi;
        pj = pi;
        pk = pi - 1;
        while (pj > pl && DTYPE_LT(vp, *pk)) {
            *pj-- = *pk--;
        }
        *pj = vp;
    }
}


void select_sort(dtype *a, int n) {
    int i, j, k;
    for (i=0; i<n-1; ++i) {
        k = i;
        for (j=i+1; j<n; ++j) {
            if (a[j] < a[k]) k = j;
        }
        DTYPE_SWAP(a[i], a[k]);
    }
}


void numpy_quicksort(dtype *start, int num)
{
    dtype vp;
    dtype *pl = start;
    dtype *pr = start + num - 1;
    dtype *stack[PYA_QS_STACK];
    dtype **sptr = stack;
    dtype *pm, *pi, *pj, *pk;

    for (;;) {
        while ((pr - pl) > SMALL_QUICKSORT) {
            /* quicksort partition */
            pm = pl + ((pr - pl) >> 1);
            if (DTYPE_LT(*pm, *pl)) DTYPE_SWAP(*pm, *pl);
            if (DTYPE_LT(*pr, *pm)) DTYPE_SWAP(*pr, *pm);
            if (DTYPE_LT(*pm, *pl)) DTYPE_SWAP(*pm, *pl);
            vp = *pm;
            pi = pl;
            pj = pr - 1;
            DTYPE_SWAP(*pm, *pj);
            for (;;) {
                do ++pi; while (DTYPE_LT(*pi, vp));
                do --pj; while (DTYPE_LT(vp, *pj));
                if (pi >= pj) {
                    break;
                }
                DTYPE_SWAP(*pi,*pj);
            }
            pk = pr - 1;
            DTYPE_SWAP(*pi, *pk);
            /* push largest partition on stack */
            if (pi - pl < pr - pi) {
                *sptr++ = pi + 1;
                *sptr++ = pr;
                pr = pi - 1;
            }
            else {
                *sptr++ = pl;
                *sptr++ = pi - 1;
                pl = pi + 1;
            }
        }

        /* insertion sort */
        for (pi = pl + 1; pi <= pr; ++pi) {
            vp = *pi;
            pj = pi;
            pk = pi - 1;
            while (pj > pl && DTYPE_LT(vp, *pk)) {
                *pj-- = *pk--;
            }
            *pj = vp;
        }
        if (sptr == stack) {
            break;
        }
        pr = *(--sptr);
        pl = *(--sptr);
    }
}


void _numpy_mergesort(dtype *pl, dtype *pr, dtype *pw)
{
    dtype vp, *pi, *pj, *pk, *pm;

    if (pr - pl > SMALL_MERGESORT) {
        /* merge sort */
        pm = pl + ((pr - pl) >> 1);
        _numpy_mergesort(pl, pm, pw);
        _numpy_mergesort(pm, pr, pw);
        for (pi = pw, pj = pl; pj < pm;) {
            *pi++ = *pj++;
        }
        pi = pw + (pm - pl);
        pj = pw;
        pk = pl;
        while (pj < pi && pm < pr) {
            if (DTYPE_LT(*pm, *pj)) {
                *pk++ = *pm++;
            }
            else {
                *pk++ = *pj++;
            }
        }
        while(pj < pi) {
            *pk++ = *pj++;
        }
    }
    else {
        /* insertion sort */
        for (pi = pl + 1; pi < pr; ++pi) {
            vp = *pi;
            pj = pi;
            pk = pi - 1;
            while (pj > pl && DTYPE_LT(vp, *pk)) {
                *pj-- = *pk--;
            }
            *pj = vp;
        }
    }
}

void numpy_mergesort(dtype *start, int num)
{
    dtype *pl, *pr, *pw;

    pl = start;
    pr = pl + num;
    pw = (dtype *) malloc((num/2) * sizeof(dtype));
    _numpy_mergesort(pl, pr, pw);
    free(pw);
}


void quick_sort (int *a, int n) {
    if (n < 2)
        return;
    int p = a[n / 2];
    int *l = a;
    int *r = a + n - 1;
    while (l <= r) {
        if (*l < p) {
            l++;
            continue;
        }
        if (*r > p) {
            r--;
            continue; // we need to check the condition (l <= r) every time we change the value of l or r
        }
        int t = *l;
        *l++ = *r;
        *r-- = t;
    }
    quick_sort(a, r - a + 1);
    quick_sort(l, a + n - l);
}


void _merge_sort(dtype *a, int n, dtype *scratch) {
    int i, m=n/2;
    dtype *l=a, *r=a+m;
    
    if (!m) return;

    _merge_sort(l, m, scratch);
    _merge_sort(r, n - m, scratch + m);

    for (i=0; i<n; ++i) {
        scratch[i] = a[i];
    }
    l = scratch;
    r = l + m;
    while (l < scratch + m && r < scratch + n) {
        *a++ = DTYPE_LT(*l, *r) ? *l++ : *r++;
    }
    if (r == scratch + n) {
        while (l < scratch + m) {
            *a++ = *l++;
        } 
    } else {
        while (r < scratch + n) {
            *a++ = *r++;
        }
    }
}

void merge_sort(dtype *a, int n) {
    dtype *scratch = malloc(n*sizeof(dtype));
    _merge_sort(a, n, scratch);
    free(scratch);
}


typedef struct {
    int x;
    int i;
} Item;

Item *a, v, p, t;

#define exch(a, b) do {t=a; a=b; b=t;} while(0)

void print_array(Item *a, int l, int r) {
    int i;
    for (i=l; i<=r; ++i) {
        printf("%d ", a[i].i);
    }
    printf("\n");
}


void quick_sort_arg(int i, int n) {
    if (n < 2) return;

    //print_array(a, i, i + n - 1);

    p = a[i + n / 2];
    int l = i, r = i + n - 1;
    while (l <= r) {
        if (a[l].x < p.x) {
            ++l;
            continue;
        }
        if (a[r].x > p.x) {
            --r;
            continue;
        }
        t = a[l];
        a[l++] = a[r];
        a[r--] = t;
    }
    quick_sort_arg(i, r - i + 1);
    quick_sort_arg(l, i + n - l);
}


int compare(Item a, Item b) {
    return a.x < b.x;
    //return (a.x < b.x) || (a.x == b.x) && (a.i < b.i);
}

void quick_sort_arg2(int l, int r) {
    if (r - l < 1) return;
    
    //print_array(a, l, r);   
    int ll=l, rr=r, mm=l + (r + 1 - l) / 2;
    p = a[mm];
    while (ll <= rr) {
        if (compare(a[ll], p) || (a[ll].x == p.x) && (ll < mm)) {
            ++ll;
            continue;
        }
        if (compare(p, a[rr]) || (p.x == a[rr].x) && (mm < rr)) {
            --rr;
            continue;
        }
        t = a[ll];
        a[ll++] = a[rr];
        a[rr--] = t;
    }
    quick_sort_arg2(l, rr);
    quick_sort_arg2(ll, r);
}



void quicksort_3way(int l, int r) { 
    int i = l-1, j = r, p = l-1, q = r; 
    Item v = a[r];
    if (r <= l) return; 
    for (;;) {
        while (a[++i].x < v.x);
        while (v.x < a[--j].x) if (j == l) break; 
        if (i >= j) break;
        exch(a[i], a[j]);
    }
    exch(a[i], a[r]); j = i-1; i = i+1;
    quicksort_3way(l, j);
    quicksort_3way(i, r); 
}


#define SIZE 1000000 /* 1 MB */

static const unsigned char masks[8] = {1, 2, 4, 8, 16, 32, 64, 128};

unsigned int
bitmap_sort(unsigned int *a, int n, unsigned int *res) {
    unsigned char bitarray[SIZE] = {0}, bitelem;
    unsigned int x, *t, i, j, k, MAX=8*SIZE;

    t = a;
    for (i=n; i--; ) {
        x = *t++;
        j = x >> 3;
        k = x - (j << 3);
        bitarray[j] |= masks[k];
    }

    if (!res) res = a;
    t = res;
    for (j=0; j<SIZE; ++j) {
        bitelem = bitarray[j];
        for (k=0; k<8; ++k) {
            if (bitelem & masks[k]) {
                *t++ = (j << 3) + k;
            }
        }
    }
    return t - res;
}


/*
int main () {
    //Item array[] = {{2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, 
    //                {2, 7}, {2, 8}, {2, 9}};

    Item array[] = {{4, 0}, {65, 1}, {2, 2}, {-31, 3}, {0, 4}, {99, 5}, {2, 6}, 
                    {83, 7}, {782, 8}, {1, 9}};
    int n = sizeof array / sizeof array[0];
    a = array;
    //quick_sort(array, n);
    //quick_sort_arg(0, n);
    //quick_sort_arg2(0, n - 1); 
    quicksort_3way(0, n - 1);
    int i;
    for (i=0; i<n; ++i) printf("%d %d\n", array[i].x, array[i].i);
    return 0;
}
*/

