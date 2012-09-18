#include "pqueue2.h"

#include "pqueue.h"

#define PQUEUE2

#ifndef PQUEUE2

void *pqueue2_init(size_t n, void *cmppri, 
        void *getpri, void *setpri, void *getpos, void *setpos) 
{
    pqueue2_t *Q = malloc(sizeof(pqueue2_t));
    Q->cmppri = cmppri;
    Q->q = pqueue_init(n, cmppri, getpri, setpri, getpos, setpos);
    return Q;
}

void pqueue2_free(pqueue2_t *Q) {
    pqueue_free(Q->q);
    free(Q);
}

void *pqueue2_replace_head(pqueue2_t *Q, void *d) {
    return pqueue_replace_head(Q->q, d);
}

void pqueue2_replace_with_higher(pqueue2_t *Q, void *n, void *d) {
    pqueue_replace_with_higher(Q->q, n, d);
}

int pqueue2_insert(pqueue2_t *Q, void *d) {
    return pqueue_insert(Q->q, d);
}

void pqueue2_update(pqueue2_t *Q, int cmp, void *d) {
    pqueue_update(Q->q, cmp, d);
}

void *pqueue2_pop(pqueue2_t *Q) {
    return pqueue_pop(Q->q);
}

int pqueue2_remove(pqueue2_t *Q, void *d) {
    return pqueue_remove(Q->q, d);
}

void *pqueue2_peek(pqueue2_t *Q) {
    return pqueue_peek(Q->q);
}

void pqueue2_print(pqueue2_t *Q, void *out, void *print) {
    pqueue_print(Q->q, out, print);
}


#else

typedef struct {
    int n_head;
    int n_bulk;
    void *head;
    void *bulk;
} dbl_queue_t;

void *pqueue2_init(size_t n_bulk, void *cmppri, 
        void *getpri, void *setpri, void *getpos, void *setpos) 
{
    static const int n_head = 16;
    pqueue2_t *Q = malloc(sizeof(pqueue2_t));
    Q->cmppri = cmppri;

    dbl_queue_t *q = malloc(sizeof(dbl_queue_t));
    q->n_head = n_head;
    q->n_bulk = n_bulk;
    q->head = pqueue_init(n_head, cmppri, getpri, setpri, getpos, setpos);
    q->bulk = pqueue_init(n_bulk, cmppri, getpri, setpri, getpos, setpos);
    Q->q = q;
    return Q;
}

void pqueue2_free(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    pqueue_free(q->bulk);
    pqueue_free(q->head);
    free(q);
    free(Q);
}

void *pqueue2_replace_head(pqueue2_t *Q, void *d) {
    //dbl_queue_t *q = Q->q;
    //return pqueue_replace_head(q->bulk, d);
    void *n = pqueue2_pop(Q);
    pqueue2_insert(Q, d);
    return n;
}

void pqueue2_replace_with_higher(pqueue2_t *Q, void *n, void *d) {
    //dbl_queue_t *q = Q->q;
    //pqueue_replace_with_higher(q->bulk, n, d);
    pqueue2_remove(Q, n);
    pqueue2_insert(Q, d);
}

int pqueue2_insert(pqueue2_t *Q, void *d) {
    dbl_queue_t *q = Q->q;
    return pqueue_insert(q->bulk, d);
}

void pqueue2_update(pqueue2_t *Q, int cmp, void *d) {
    dbl_queue_t *q = Q->q;
    pqueue_update(q->bulk, cmp, d);
}

void *pqueue2_pop(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    return pqueue_pop(q->bulk);
}

int pqueue2_remove(pqueue2_t *Q, void *d) {
    dbl_queue_t *q = Q->q;
    return pqueue_remove(q->bulk, d);
}

void *pqueue2_peek(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    return pqueue_peek(q->bulk);
}

void pqueue2_print(pqueue2_t *Q, void *out, void *print) {
    dbl_queue_t *q = Q->q;
    pqueue_print(q->bulk, out, print);
}


#endif




