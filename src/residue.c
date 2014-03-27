#include <stdlib.h>
#include "residue.h"

struct residue *residue_alloc(const struct AA *aa){
    struct residue *r = malloc(sizeof(struct residue));
    if(!r)
        return NULL;
    residue_init(r, aa);
    return r;
}

void residue_init(struct residue *r, const struct AA *aa){
    r->aa = aa;
    vector_zero(&r->position);
    vector_zero(&r->velocity);
    vector_zero(&r->force);
}

void residue_free(struct residue *r){
    free(r);
}
