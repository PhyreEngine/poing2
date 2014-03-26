#include <stdlib.h>
#include "residue.h"

struct residue *residue_alloc(const struct AA *aa){
    struct residue *r = malloc(sizeof(struct residue));
    if(!r)
        return NULL;
    r->aa = aa;
    vector_zero(&r->position);
    vector_zero(&r->velocity);
    vector_zero(&r->force);
    return r;
}

void residue_free(struct residue *r){
    free(r);
}
