#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "vector.h"

struct residue {
    const struct AA *aa;
    struct vector position, velocity, force;
};

struct residue *residue_alloc(const struct AA *aa);
void residue_free(struct residue *r);
extern struct AA * AA_lookup (register const char *str, register unsigned int len);

#endif /* RESIDUE_H_ */
