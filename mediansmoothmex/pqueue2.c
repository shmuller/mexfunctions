#include "pqueue2.h"

#include "pqueue.h"

#ifndef PQUEUE2

void *pqueue2_init(size_t n, void *cmppri, 
        void *getpri, void *setpri, void *getpos, void *setpos) 
{
    pqueue2_t *q = malloc(sizeof(pqueue2_t));
    q->cmppri = cmppri;
    q->q = pqueue_init(n, cmppri, getpri, setpri, getpos, setpos);
    return q;
}

void pqueue2_free(pqueue2_t *q) {
    pqueue_free(q->q);
    free(q);
}

void *pqueue2_replace_head(pqueue2_t *q, void *d) {
    return pqueue_replace_head(q->q, d);
}

void pqueue2_replace_with_higher(pqueue2_t *q, void *n, void *d) {
    pqueue_replace_with_higher(q->q, n, d);
}

int pqueue2_insert(pqueue2_t *q, void *d) {
    return pqueue_insert(q->q, d);
}

void pqueue2_update(pqueue2_t *q, int cmp, void *d) {
    pqueue_update(q->q, cmp, d);
}

void *pqueue2_pop(pqueue2_t *q) {
    return pqueue_pop(q->q);
}

int pqueue2_remove(pqueue2_t *q, void *d) {
    return pqueue_remove(q->q, d);
}

void *pqueue2_peek(pqueue2_t *q) {
    return pqueue_peek(q->q);
}

void pqueue2_print(pqueue2_t *q, void *out, void *print) {
    pqueue_print(q->q, out, print);
}


#else

#endif




