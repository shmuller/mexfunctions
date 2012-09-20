#include "pqueue2.h"

#include "pqueue.h"

#include <stdio.h>

typedef struct {
    int n_head;
    int n_bulk;
    void *head;
    void *bulk;
} dbl_queue_t;


int insert_head(dbl_queue_t *q, void *d) {
    ++count.insert_head;
    ((node_t*)d)->id2 = 0;
    return pqueue_insert(q->head, d);
}

int insert_bulk(dbl_queue_t *q, void *d) {
    ++count.insert_bulk;
    ((node_t*)d)->id2 = 1;
    return pqueue_insert(q->bulk, d);
}

int remove_head(dbl_queue_t *q, void *d) {
    ++count.remove_head;
    return pqueue_remove(q->head, d);
}

int remove_bulk(dbl_queue_t *q, void *d) {
    ++count.remove_bulk;
    return pqueue_remove(q->bulk, d);
}


void *pqueue2_init(unsigned char id, size_t n_bulk, void *cmppri, 
    void *getpri, void *setpri, void *getpos, void *setpos) 
{
    static const int n_head = 4;
    pqueue2_t *Q = malloc(sizeof(pqueue2_t));
    Q->cmppri = cmppri;
    Q->getpri = getpri;
    Q->id = id;

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
    ++count.replace_head;
    void *n = pqueue2_pop(Q);
    pqueue2_insert(Q, d);
    return n;
}

void pqueue2_replace_with_higher(pqueue2_t *Q, void *n, void *d) {
    //dbl_queue_t *q = Q->q;
    //pqueue_replace_with_higher(q->bulk, n, d);
    ++count.replace_with_higher;
    pqueue2_remove(Q, n);
    pqueue2_insert(Q, d);
}

int pqueue2_insert(pqueue2_t *Q, void *d) {
    dbl_queue_t *q = Q->q;
    void *b, *h;

    b = pqueue_peek(q->bulk);
    if (!b || Q->cmppri(Q->getpri(b), Q->getpri(d))) {
        // new element ranks higher than first element of bulk
        if (pqueue_size(q->head) < q->n_head) {
            // head is not full
            insert_head(q, d);
        } else {
            // head is full, need to rebalance
            h = pqueue_peek_last(q->head);
            if (h && Q->cmppri(Q->getpri(h), Q->getpri(d))) {
                // new element ranks higher than last element of head
                //remove_head(q, h);
                //insert_head(q, d);
                ((node_t*)d)->id2 = 0;
                pqueue_replace_with_higher(q->head, h, d);
                insert_bulk(q, h);
            } else {
                // new element belongs to top of bulk
                insert_bulk(q, d);
            }
        }        
    } else {
        // new element ranks lower than first element of the bulk
        insert_bulk(q, d);
    }

    return 0;
}

void pqueue2_update(pqueue2_t *Q, int cmp, void *d) {
    dbl_queue_t *q = Q->q;
    int is_wrong;
    if (((node_t*)d)->id2 == 0) {
        pqueue_update(q->head, cmp, d);
        is_wrong = (d == pqueue_peek_last(q->head));
    } else {
        pqueue_update(q->bulk, cmp, d);
        is_wrong = (d == pqueue_peek(q->bulk));
    }

    if (is_wrong) {
        ++count.update_is_wrong;
        pqueue2_remove(Q, d);
        pqueue2_insert(Q, d);
    }
}

void *pqueue2_pop(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    void *head = pqueue_pop(q->head);
    if (!head)
        head = pqueue_pop(q->bulk);
    return head;
}

int pqueue2_remove(pqueue2_t *Q, void *d) {
    dbl_queue_t *q = Q->q;
    if (((node_t*)d)->id2 == 0)
        return remove_head(q, d);
    else
        return remove_bulk(q, d);
}

void *pqueue2_peek(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    void *head = pqueue_peek(q->head);
    if (!head)
        head = pqueue_peek(q->bulk);
    return head;
}

void *pqueue2_peek_last(pqueue2_t *Q) {
    dbl_queue_t *q = Q->q;
    void *last = pqueue_peek_last(q->bulk);
    if (!last)
        last = pqueue_peek_last(q->head);
    return last;
}

void pqueue2_print(pqueue2_t *Q, void *out, void *print) {
    dbl_queue_t *q = Q->q;
    pqueue_print(q->head, out, print);
    printf("--");
    pqueue_print(q->bulk, out, print);
}

void pqueue2_dump(pqueue2_t *Q, void *out, void *print) {
    dbl_queue_t *q = Q->q;
    pqueue_dump(q->head, out, print);
    pqueue_dump(q->bulk, out, print);
}


