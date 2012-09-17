/* Sliding median filter implementation using pqueues.
 *
 * S. H. Muller, 2012/09/13
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mediansmooth_pqueue.h"

#include "pqueue.h"

// definitions for pqueue
typedef struct {
	void *pri;
	size_t pos;
    pqueue_t *q;
} node_t;

static void *get_pri(void *a) {
	return ((node_t*)a)->pri;
}

static void set_pri(void *a, void *pri) {
    ((node_t*)a)->pri = pri;
}

static size_t get_pos(void *a) {
	return ((node_t*)a)->pos;
}

static void set_pos(void *a, size_t pos) {
	((node_t*)a)->pos = pos;
}


// operations on pqueues
void pqueue_insert_wrap(pqueue_t *q, node_t *n) {
    n->q = q;
    pqueue_insert(q, n);
}

void *pqueue_replace_head_wrap(pqueue_t *q, node_t *n) {
    n->q = q;
    return pqueue_replace_head(q, n);
}

void pqueue_replace_with_higher_wrap(pqueue_t *q, node_t *o, node_t *n) {
    n->q = q;
    pqueue_replace_with_higher(q, o, n);
}

void insert(pqueue_t *dest, pqueue_t *src, node_t *n, void *med_pri) {
    int cmp = src->cmppri(n->pri, med_pri);
    node_t *o = (cmp) ? pqueue_replace_head_wrap(src, n) : n;
    pqueue_insert_wrap(dest, o);
}

void delete(pqueue_t *q_long, pqueue_t *q_short, node_t *n) {
    if (n->q == q_long)
        pqueue_remove(q_long, n);
    else
        pqueue_replace_with_higher_wrap(q_short, n, pqueue_pop(q_long));
}

void pqueue_printfun(FILE *out, void *a) {
    ((printfun_t*)out)(get_pri(a));
    printf(", ");
}

void print_queues_median(pqueue_t *L, pqueue_t *R, void *med, printfun_t *print) {
    pqueue_print(L, (FILE*) print, &pqueue_printfun);
    printf("|"); print(med); printf("|, ");
    pqueue_print(R, (FILE*) print, &pqueue_printfun);
    printf("\n");
}


// median status and update operations
typedef struct {
    int balance;
    node_t *median;
    pqueue_t *L;
    pqueue_t *R;
} status;

void status_init(status *S, int m, void *lt, void *gt) {
    S->balance = 0;
    S->median = NULL;
    S->L = pqueue_init(m, lt, get_pri, set_pri, get_pos, set_pos);
	S->R = pqueue_init(m, gt, get_pri, set_pri, get_pos, set_pos);    
}

void status_free(status *S) {
    pqueue_free(S->L);
    pqueue_free(S->R);
}

void *nodes_init(int W, int BYTES) {
    int i;
    void *nodes = malloc(W*BYTES);
    node_t *n = nodes;
    for (i=0; i < W; ++i) {
        n = nodes + i*BYTES;
        n->pri = n + 1;
    }
    return nodes;
}

void *add(status *S, node_t *n, void *new_pri, int bytes) {
    memcpy(n->pri, new_pri, bytes);
    
    if (!S->median)
        S->median = n;

    switch (S->balance) {
        case 0:
            if (S->R->cmppri(n->pri, S->median->pri)) {
                pqueue_insert_wrap(S->R, n);
                S->balance = +1;
                S->median = pqueue_peek(S->R);
            } else {
                pqueue_insert_wrap(S->L, n);
                S->balance = -1;
                S->median = pqueue_peek(S->L);
            }
            break;
        case +1:
            insert(S->L, S->R, n, S->median->pri);
            S->balance = 0;
            S->median = pqueue_peek(S->L);
            break;
        case -1:
            insert(S->R, S->L, n, S->median->pri);
            S->balance = 0;
            S->median = pqueue_peek(S->L);
            break;
    }
    return S->median->pri;
}

void *del(status *S, node_t *n) {
    switch (S->balance) {
        case 0:
            if (n->q == S->L) {
                pqueue_remove(S->L, n);
                S->balance = +1;
                S->median = pqueue_peek(S->R);
            } else {
                pqueue_remove(S->R, n);
                S->balance = -1;
                S->median = pqueue_peek(S->L);
            }
            break;
        case +1:
            delete(S->R, S->L, n);
            S->balance = 0;
            S->median = pqueue_peek(S->L);
            break;
        case -1:
            delete(S->L, S->R, n);
            S->balance = 0;
            S->median = pqueue_peek(S->L);
            break;
    }
    return S->median->pri;
}

void *rep(status *S, node_t *n, void *new_pri, int bytes) {
    int cmp;
    if (n->q->cmppri(new_pri, S->median->pri)) {
        // shortcut: new and old elements belong to the same queue
        cmp = n->q->cmppri(n->pri, new_pri);
        memcpy(n->pri, new_pri, bytes);
        pqueue_update(n->q, cmp, n);
        // median comes from the same queue, but may have shifted
        S->median = pqueue_peek(S->median->q);
    } else {
        del(S, n);
        add(S, n, new_pri, bytes);
    }
    return S->median->pri;
}


// actual median filter implementations for different boundary treatments

#define NEXT_NODE (nodes + (j % W)*BYTES)

void median_filt_pqueue_bdry_0(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x, *y;
    int i, j=0, W = 2*w + 1, BYTES = sizeof(node_t) + bytes;

    void *nodes = nodes_init(W, BYTES);
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    x = y = X;
    for (i=0; i < w; ++i, ++j, x += bytes) {
        med = add(&S, NEXT_NODE, x, bytes);
    }
    for (i=0; i < w+1; ++i, ++j, x += bytes, y += bytes) {
        med = add(&S, NEXT_NODE, x, bytes);
        memcpy(y, med, bytes);
    }
    for (; i < N-w; ++i, ++j, x += bytes, y += bytes) {
        med = rep(&S, NEXT_NODE, x, bytes);
        memcpy(y, med, bytes);
    }
    for (; i < N; ++i, ++j, y += bytes) {
        med = del(&S, NEXT_NODE);
        memcpy(y, med, bytes);
    }

    status_free(&S);
    free(nodes);
}


void median_filt_pqueue_bdry_1(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x, *y;
    int i=1, j=0, W = 2*w + 1, BYTES = sizeof(node_t) + bytes;
    
    void *nodes = nodes_init(W, BYTES);
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    x = y = X;
    med = add(&S, NEXT_NODE, x, bytes); ++j; x += bytes;
    memcpy(y, med, bytes); y += bytes;

    for (; i < w+1; ++i) {
        med = add(&S, NEXT_NODE, x, bytes); ++j; x += bytes;
        med = add(&S, NEXT_NODE, x, bytes); ++j; x += bytes;
        memcpy(y, med, bytes); y += bytes;
    }
    for (; i < N-w; ++i) {
        med = rep(&S, NEXT_NODE, x, bytes); ++j; x += bytes;
        memcpy(y, med, bytes); y += bytes;
    }
    for (; i < N; ++i) {
        med = del(&S, NEXT_NODE); ++j;
        med = del(&S, NEXT_NODE); ++j;
        memcpy(y, med, bytes); y += bytes;
    }

    status_free(&S);
    free(nodes);
}


void median_filt_pqueue_bdry_2(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x = X;
    int i, BYTES = sizeof(node_t) + bytes;
    
    void *nodes = nodes_init(w, BYTES);
    status S;
    status_init(&S, (w+1)/2, fun->lt, fun->gt);

    for (i=0; i<w; ++i, x += bytes) {
        med = add(&S, nodes + i*BYTES, x, bytes);
        memcpy(x, med, bytes);
        //print_queues_median(S.L, S.R, med, fun->print);
        //printf("L=%d, R=%d\n", (int) pqueue_size(S.L), (int) pqueue_size(S.R));
    }

    for(; i<N; ++i, x += bytes) {
        med = rep(&S, nodes + (i % w)*BYTES, x, bytes);
        memcpy(x, med, bytes);
        //print_queues_median(S.L, S.R, med, fun->print);
        //printf("L=%d, R=%d\n", (int) pqueue_size(S.L), (int) pqueue_size(S.R));
    }

    status_free(&S);
    free(nodes);
}


double mean(void *y, int N, int bytes, to_double_t *to_double) {
    double m = 0.;
    int i;
    for (i=N; i--; ) {
        m += to_double(y);
        y += bytes;
    }
    return m / N;
}

double slope(void *y, int N, int bytes, to_double_t *to_double) {
    double k = 0.;
    int i, M = N-1;
    for (i=-M; i<=M; i+=2) {
        k += i * to_double(y);
        y += bytes;
    }
    return k * 6. / M / N / (N+1);
}

void linfit(double *kd, void *x, int n, int bytes, to_double_t *to_double) {
    double m = mean(x, n, bytes, to_double);
    kd[0] = slope(x, n, bytes, to_double);
    kd[1] = m - kd[0]*(n-1)/2.;
}


void median_filt_pqueue_bdry_3(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x, *y, *buf = malloc(bytes);
    int i, j=0, W = 2*w + 1, BYTES = sizeof(node_t) + bytes;
    double kd0[2], kd1[2];
    
    linfit(kd0, X, w, bytes, fun->to_double);
    linfit(kd1, X + (N-w)*bytes, w, bytes, fun->to_double);

    void *nodes = nodes_init(W, BYTES);
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    x = X;
    for (i=0; i < w; ++i, ++j, x += bytes) {
        med = add(&S, NEXT_NODE, x, bytes);
        fun->from_double(x, kd0[0]*i + kd0[1]);
    }
    y = x;
    for (i=0; i < w+1; ++i, ++j, x += bytes) {
        med = add(&S, NEXT_NODE, x, bytes);
    }
    memcpy(y, med, bytes);
    y += bytes;

    for (; i < N-w; ++i, ++j, x += bytes, y += bytes) {
        med = rep(&S, NEXT_NODE, x, bytes);
        memcpy(y, med, bytes);
    }
    for (i=0; i < w; ++i, y += bytes) {
        fun->from_double(y, kd1[0]*i + kd1[1]);
    }

    status_free(&S);
    free(nodes);
    free(buf);
}


void median_filt_pqueue_bdry_4(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x, *y, *buf = malloc(bytes);
    int i, j=0, W = 2*w + 1, BYTES = sizeof(node_t) + bytes;
    double kd0[2], kd1[2];
    
    linfit(kd0, X, w + 1, bytes, fun->to_double);
    linfit(kd1, X + (N-w-1)*bytes, w + 1, bytes, fun->to_double);

    void *nodes = nodes_init(W, BYTES);
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    x = y = X;
    for (i=-w; i < 0; ++i, ++j) {
        fun->from_double(buf, kd0[0]*i + kd0[1]);
        med = add(&S, NEXT_NODE, buf, bytes);
    }
    for (; i < w+1; ++i, ++j, x += bytes) {
        med = add(&S, NEXT_NODE, x, bytes);
    }
    memcpy(y, med, bytes);
    y += bytes;

    for (; i < N; ++i, ++j, x += bytes, y += bytes) {
        med = rep(&S, NEXT_NODE, x, bytes);
        memcpy(y, med, bytes);
    }
    for (i=0; i < w; ++i, ++j, y += bytes) {
        fun->from_double(buf, kd1[0]*(i+w+1) + kd1[1]);
        med = rep(&S, NEXT_NODE, buf, bytes);
        memcpy(y, med, bytes);
    }
    
    status_free(&S);
    free(nodes);
    free(buf);
}


void median_filt_pqueue(void *X, int N, int w, int *ind, int bdry, int bytes, fun_t *fun) {
    switch (bdry) {
        case 0:
            median_filt_pqueue_bdry_0(X, N, w, ind, bytes, fun);
            break;
        case 1:
            median_filt_pqueue_bdry_1(X, N, w, ind, bytes, fun);
            break;
        case 2:
            median_filt_pqueue_bdry_2(X, N, w, ind, bytes, fun);
            break;
        case 3:
            median_filt_pqueue_bdry_3(X, N, w, ind, bytes, fun);
            break;
        case 4:
            median_filt_pqueue_bdry_4(X, N, w, ind, bytes, fun);
            break;

    }
}


