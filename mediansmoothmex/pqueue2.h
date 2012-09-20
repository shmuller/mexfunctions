
#include <stdlib.h> /* size_t */

typedef struct {
	void *pri;
	size_t pos;
    unsigned char id;
    unsigned char id2;
} node_t;

typedef int (pqueue2_cmp_pri_f)(void *, void *);
typedef void* (pqueue2_get_pri_f)(void *a);

typedef struct {
    pqueue2_cmp_pri_f *cmppri;
    pqueue2_get_pri_f *getpri;
    unsigned char id;
    void *q;
} pqueue2_t;

void *pqueue2_init(unsigned char id, size_t n, void *cmppri, 
    void *getpri, void *setpri, void *getpos, void *setpos);

void pqueue2_free(pqueue2_t *q);
void *pqueue2_replace_head(pqueue2_t *q, void *d);
void pqueue2_replace_with_higher(pqueue2_t *q, void *n, void *d);
int pqueue2_insert(pqueue2_t *q, void *d);
void pqueue2_update(pqueue2_t *q, int cmp, void *d);
void *pqueue2_pop(pqueue2_t *q);
int pqueue2_remove(pqueue2_t *q, void *d);
void *pqueue2_peek(pqueue2_t *q);
void *pqueue2_peek_last(pqueue2_t *q);
void pqueue2_print(pqueue2_t *q, void *out, void *print);
void pqueue2_dump(pqueue2_t *q, void *out, void *print);

// debugging
typedef struct {
    int add_0;
    int add_1;
    int del_0;
    int del_1;
    int update;
    int update_is_wrong;
    int insert_head;
    int insert_bulk;
    int remove_head;
    int remove_bulk;
    int replace_head;
    int replace_with_higher;
} count_t;

count_t count;

void pqueue2_stats_reset();
void pqueue2_stats_print();

