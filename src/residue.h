#ifndef RESIDUE_H_
#define RESIDUE_H_

#include <stdbool.h>
#include "vector.h"

struct AA {
    const char *oneletter;
    const char *threeletter;
};

struct residue {
    const struct AA *aa;
    int id;
    bool synthesised;
    struct vector position, velocity, force;
};

struct residue *residue_alloc(const struct AA *aa, int id);
void residue_init(struct residue *r, const struct AA *aa, int id);
void residue_free(struct residue *r);
extern struct AA * AA_lookup (register const char *str, register unsigned int len);

#endif /* RESIDUE_H_ */
