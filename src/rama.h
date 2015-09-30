#ifndef RAMA_H_
#define RAMA_H_

#include <stdbool.h>
#include "residue.h"
#include "torsion_spring.h"

enum rama_constraint_type {
    GENERAL,
    ALPHA,
    BETA,
    PROLINE,
    PRE_PROLINE,
    GLYCINE,
    ALANINE,
    UNKNOWN_RAMA
};

struct rama_constraint {
    struct torsion_spring *phi, *psi;
    enum rama_constraint_type type;
    bool enabled;
};

int rama_read_closest(const char *file, enum rama_constraint_type type);
int rama_get_closest(struct rama_constraint *rama);
void rama_free_data();
enum rama_constraint_type rama_parse_type(const char *type);
void rama_init(struct rama_constraint *rama,
        struct residue *residue,
        struct residue *next_residue,
        const char *type);
void rama_random_init(struct rama_constraint *rama);

#endif //RAMA_H_
