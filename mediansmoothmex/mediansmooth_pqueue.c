/* Sliding median filter implementation using pqueues.
 *
 * S. H. Muller, 2012/09/13
 */

#include <stdio.h>
#include <stdlib.h>

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
    // never used!
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

void *add(status *S, node_t *n, void *new_pri) {
    n->pri = new_pri;

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

void *rep(status *S, node_t *n, void *new_pri) {
    void *old_pri;
    if (n->q->cmppri(new_pri, S->median->pri)) {
        // shortcut: new and old elements belong to the same queue
        old_pri = n->pri;
        n->pri = new_pri;
        pqueue_change_priority2(n->q, old_pri, n);
        // median comes from the same queue, but may have shifted
        S->median = pqueue_peek(S->median->q);
    } else {
        del(S, n);
        add(S, n, new_pri);
    }
    return S->median->pri;
}


// actual median filter implementations for different boundary treatments

#define NEXT_NODE (nodes + (j++ % W))

void median_filt_pqueue_bdry_0(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x = X;
    int i, j=0, W = 2*w + 1;
    
    node_t *nodes = malloc(W*sizeof(node_t));
	
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    for (i=0; i < w; ++i) {
        med = add(&S, NEXT_NODE, x); x += bytes;
    }
    for (i=0; i < w+1; ++i) {
        med = add(&S, NEXT_NODE, x); x += bytes;
        *ind++ = (med - X) / bytes;
    }
    for (; i < N-w; ++i) {
        med = rep(&S, NEXT_NODE, x); x += bytes;
        *ind++ = (med - X) / bytes;
    }
    for (; i < N; ++i) {
        med = del(&S, NEXT_NODE);
        *ind++ = (med - X) / bytes;
    }

    status_free(&S);
    free(nodes);
}


void median_filt_pqueue_bdry_1(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x = X;
    int i=1, j=0, W = 2*w + 1;
    
    node_t *nodes = malloc(W*sizeof(node_t));
	
    status S;
    status_init(&S, w + 1, fun->lt, fun->gt);

    med = add(&S, NEXT_NODE, x); x += bytes;
    *ind++ = (med - X) / bytes;

    for (; i < w+1; ++i) {
        med = add(&S, NEXT_NODE, x); x += bytes;
        med = add(&S, NEXT_NODE, x); x += bytes;
        *ind++ = (med - X) / bytes;
    }
    for (; i < N-w; ++i) {
        med = rep(&S, NEXT_NODE, x); x += bytes;
        *ind++ = (med - X) / bytes;
    }
    for (; i < N; ++i) {
        med = del(&S, NEXT_NODE);
        med = del(&S, NEXT_NODE);
        *ind++ = (med - X) / bytes;
    }

    status_free(&S);
    free(nodes);
}


void median_filt_pqueue_bdry_2(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    void *med, *x = X;
    int i;
    
    node_t *nodes = malloc(w*sizeof(node_t));
	
    status S;
    status_init(&S, (w+1)/2, fun->lt, fun->gt);

    for (i=0; i<w; ++i) {
        med = add(&S, nodes + i, x); x += bytes;
        ind[i] = (med - X) / bytes;
        print_queues_median(S.L, S.R, med, fun->print);
        //printf("L=%d, R=%d\n", (int) pqueue_size(S.L), (int) pqueue_size(S.R));
    }

    for(; i<N; ++i) {
        med = rep(&S, nodes + (i % w), x); x += bytes;
        ind[i] = (med - X) / bytes;
        print_queues_median(S.L, S.R, med, fun->print);
        //printf("L=%d, R=%d\n", (int) pqueue_size(S.L), (int) pqueue_size(S.R));
    }

    status_free(&S);
    free(nodes);
}


double mean(void *y, int N, int bytes, as_double_t *as_double) {
    double m = 0.;
    int i;
    for (i=N; i--; ) {
        m += as_double(y);
        y += bytes;
    }
    return m / N;
}

double slope(void *y, int N, int bytes, as_double_t *as_double) {
    double k = 0.;
    int i, M = N-1;
    for (i=-M; i<=M; i+=2) {
        k += i * as_double(y);
        y += bytes;
    }
    return k * 6. / M / N / (N+1);
}

void median_filt_pqueue_bdry_3(void *X, int N, int w, int *ind, int bytes, fun_t *fun) {
    double m = mean(X, N, bytes, fun->as_double);
    double k = slope(X, N, bytes, fun->as_double);
    int i;

    for (i=0; i<N; ++i) {
        ind[i] = k*(i - (N-1)/2.) + m;
    }

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
    }
}


