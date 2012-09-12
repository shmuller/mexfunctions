#include <stdio.h>
#include <stdlib.h>

#include "pqueue.h"

typedef double dtype;

typedef struct {
	void *pri;
	size_t pos;
    pqueue_t *q;
} node_t;

void pqueue_insert_wrap(pqueue_t *q, node_t *n) {
    n->q = q;
    pqueue_insert(q, n);
}

void *pqueue_replace_head_wrap(pqueue_t *q, node_t *n) {
    n->q = q;
    return pqueue_replace_head(q, n);
}


static int compare(void *next, void *curr) {
    dtype a = *(dtype*)next;
    dtype b = *(dtype*)curr;
    return (a>b) - (a<b);
}

static int lt(void *next, void *curr) {
	return *(dtype*)next < *(dtype*)curr;
}

static int gt(void *next, void *curr) {
	return *(dtype*)next > *(dtype*)curr;
}

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

void printfun(FILE *out, void *a) {
    printf("%f, ", *(dtype*)get_pri(a));
}

void print_queues_median(pqueue_t *L, pqueue_t *R, void *med) {
    pqueue_print(L, stdout, &printfun);
    printf("|%f|, ", *(dtype*)med);
    pqueue_print(R, stdout, &printfun);
    printf("\n");
}


static int balance = 0;
static node_t *median = NULL;

void insert(pqueue_t *dest, pqueue_t *src, node_t *n) {
    int cmp = src->cmppri(n->pri, median->pri);
    node_t *o = (cmp) ? pqueue_replace_head_wrap(src, n) : n;
    pqueue_insert_wrap(dest, o);
}

void *add(pqueue_t *L, pqueue_t *R, node_t *n, void *new_pri) {
    n->pri = new_pri;
    
    if (!median)
        median = n;

    switch (balance) {
        case 0:
            if (R->cmppri(n->pri, median->pri)) {
                pqueue_insert_wrap(R, n);
                balance = +1;
                median = pqueue_peek(R);
            } else {
                pqueue_insert_wrap(L, n);
                balance = -1;
                median = pqueue_peek(L);
            }
            break;
        case +1:
            insert(L, R, n);
            balance = 0;
            median = pqueue_peek(L);
            break;
        case -1:
            insert(R, L, n);
            balance = 0;
            median = pqueue_peek(L);
            break;
    }
    return median->pri;
}

void *del(pqueue_t *L, pqueue_t *R, node_t *n) {
    switch (balance) {
        case 0:
            if (n->q == L) {
                pqueue_remove(L, n);
                balance = +1;
                median = pqueue_peek(R);
            } else {
                pqueue_remove(R, n);
                balance = -1;
                median = pqueue_peek(L);
            }
            break;
        case +1:
            break;
        case -1:
            break;
    }
    return median->pri;
}

void *rep(pqueue_t *L, pqueue_t *R, node_t *n, void *new_pri) {
    if (n->q->cmppri(new_pri, median->pri)) {
        // shortcut: new and old elements belong to the same queue
        pqueue_change_priority(n->q, new_pri, n);
    } else {
        del(L, R, n);
        add(L, R, n, new_pri);
    }
    return median->pri;
}


void median_filt_pqueue(dtype *x, int N, int w, int *ind) {
	pqueue_t *L, *R;
	node_t *nodes, *n;
    dtype *med;
    int i;
    
    nodes = malloc(w*sizeof(node_t));
	L = pqueue_init(w/2+1, lt, get_pri, set_pri, get_pos, set_pos);
	R = pqueue_init(w/2+1, gt, get_pri, set_pri, get_pos, set_pos);    

    for (i=0; i<w; ++i) {
        n = &nodes[i];
        med = add(L, R, n, x + i);
        ind[i] = med - x;
        //print_queues_median(L, R, med);
    }

    for(; i<N; ++i) {
        n = &nodes[i % w];
        med = rep(L, R, n, x + i);
        ind[i] = med - x;
        //print_queues_median(L, R, med);
    }
   
	pqueue_free(L);
    pqueue_free(R);
    free(nodes);
}


