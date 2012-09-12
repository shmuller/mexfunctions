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


typedef struct {
    int balance;
    node_t *median;
} status;


void insert(pqueue_t *dest, pqueue_t *src, node_t *n, void *med_pri) {
    int cmp = src->cmppri(n->pri, med_pri);
    node_t *o = (cmp) ? pqueue_replace_head_wrap(src, n) : n;
    pqueue_insert_wrap(dest, o);
}

void *add(pqueue_t *L, pqueue_t *R, node_t *n, void *new_pri, status *stat) {
    n->pri = new_pri;

    if (!stat->median)
        stat->median = n;

    switch (stat->balance) {
        case 0:
            if (R->cmppri(n->pri, stat->median->pri)) {
                pqueue_insert_wrap(R, n);
                stat->balance = +1;
                stat->median = pqueue_peek(R);
            } else {
                pqueue_insert_wrap(L, n);
                stat->balance = -1;
                stat->median = pqueue_peek(L);
            }
            break;
        case +1:
            insert(L, R, n, stat->median->pri);
            stat->balance = 0;
            stat->median = pqueue_peek(L);
            break;
        case -1:
            insert(R, L, n, stat->median->pri);
            stat->balance = 0;
            stat->median = pqueue_peek(L);
            break;
    }
    return stat->median->pri;
}

void *del(pqueue_t *L, pqueue_t *R, node_t *n, status *stat) {
    switch (stat->balance) {
        case 0:
            if (n->q == L) {
                pqueue_remove(L, n);
                stat->balance = +1;
                stat->median = pqueue_peek(R);
            } else {
                pqueue_remove(R, n);
                stat->balance = -1;
                stat->median = pqueue_peek(L);
            }
            break;
        case +1:
            break;
        case -1:
            break;
    }
    return stat->median->pri;
}

void *rep(pqueue_t *L, pqueue_t *R, node_t *n, void *new_pri, status *stat) {
    if (n->q->cmppri(new_pri, stat->median->pri)) {
        // shortcut: new and old elements belong to the same queue
        pqueue_change_priority(n->q, new_pri, n);
        stat->median = pqueue_peek(L);
    } else {
        del(L, R, n, stat);
        add(L, R, n, new_pri, stat);
    }
    return stat->median->pri;
}


void median_filt_pqueue(dtype *x, int N, int w, int *ind) {
	pqueue_t *L, *R;
	node_t *nodes, *n;
    dtype *med;
    int i;

    status stat = {0, NULL};
    
    nodes = malloc(w*sizeof(node_t));
	L = pqueue_init(w/2+1, lt, get_pri, set_pri, get_pos, set_pos);
	R = pqueue_init(w/2+1, gt, get_pri, set_pri, get_pos, set_pos);    

    for (i=0; i<w; ++i) {
        n = &nodes[i];
        med = add(L, R, n, x + i, &stat);
        ind[i] = med - x;
        //print_queues_median(L, R, med);
        //printf("L=%d, R=%d\n", (int) pqueue_size(L), (int) pqueue_size(R));
    }

    for(; i<N; ++i) {
        n = &nodes[i % w];
        med = rep(L, R, n, x + i, &stat);
        ind[i] = med - x;
        //print_queues_median(L, R, med);
        //printf("L=%d, R=%d\n", (int) pqueue_size(L), (int) pqueue_size(R));
    }

    pqueue_free(L);
    pqueue_free(R);
    free(nodes);
}


