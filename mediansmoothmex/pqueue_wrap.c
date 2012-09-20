#include "pqueue2.h"

#include "pqueue.h"

void *pqueue2_init(unsigned char id, size_t n, void *cmppri, 
    void *getpri, void *setpri, void *getpos, void *setpos) 
{
    pqueue2_t *Q = malloc(sizeof(pqueue2_t));
    Q->cmppri = cmppri;
    Q->getpri = getpri;
    Q->id = id;
    Q->q = pqueue_init(n, cmppri, getpri, setpri, getpos, setpos);
    return Q;
}

void pqueue2_free(pqueue2_t *Q) {
    pqueue_free(Q->q);
    free(Q);
}

void *pqueue2_replace_head(pqueue2_t *Q, void *d) {
    ++count.replace_head;
    return pqueue_replace_head(Q->q, d);
    //void *n = pqueue2_pop(Q);
    //pqueue2_insert(Q, d);
    //return n;
}

void pqueue2_replace_with_higher(pqueue2_t *Q, void *n, void *d) {
    ++count.replace_with_higher;
    pqueue_replace_with_higher(Q->q, n, d);
    //pqueue2_remove(Q, n);
    //pqueue2_insert(Q, d);
}

int pqueue2_insert(pqueue2_t *Q, void *d) {
    ++count.insert_head;
    return pqueue_insert(Q->q, d);
}

void pqueue2_update(pqueue2_t *Q, int cmp, void *d) {
    pqueue_update(Q->q, cmp, d);
}

void *pqueue2_pop(pqueue2_t *Q) {
    return pqueue_pop(Q->q);
}

int pqueue2_remove(pqueue2_t *Q, void *d) {
    ++count.remove_head;
    return pqueue_remove(Q->q, d);
}

void *pqueue2_peek(pqueue2_t *Q) {
    return pqueue_peek(Q->q);
}

void *pqueue2_peek_last(pqueue2_t *Q) {
    return pqueue_peek_last(Q->q);
}

void pqueue2_print(pqueue2_t *Q, void *out, void *print) {
    pqueue_print(Q->q, out, print);
}

void pqueue2_dump(pqueue2_t *Q, void *out, void *print) {
    pqueue_dump(Q->q, out, print);
}


