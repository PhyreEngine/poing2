#ifndef RAMA_H_
#define RAMA_H_

#include <stdbool.h>
#include "residue.h"
#include "torsion_spring.h"

#define DEFAULT_RAMA_CONST 0.5

struct model;

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

int rama_is_inited(enum rama_constraint_type type);
int rama_read_closest(const char *file, enum rama_constraint_type type);
int rama_get_closest(struct rama_constraint *rama);
void rama_free_data();
enum rama_constraint_type rama_parse_type(const char *type);
void rama_init(struct rama_constraint *rama,
        const struct model *m,
        size_t residue_idx,
        const char *type,
        float constant);
void rama_random_init(struct rama_constraint *rama);
int rama_is_synthesised(struct rama_constraint *rama);

#endif //RAMA_H_
