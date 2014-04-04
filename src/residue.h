#ifndef RESIDUE_H_
#define RESIDUE_H_

#include <stddef.h>
#include <stdbool.h>
#include "vector.h"

///Maximum length (including nul byte) of atom name.
#define MAX_ATOM_NAME_SZ 5

///CA-CA bond length
#define CA_CA_LEN 3.8

///Radius of steric sphere for CA atoms
#define CA_STERIC_RADIUS 1.0

struct AA {
    const char *oneletter;
    const char *threeletter;
    double mass;
    bool has_sidechain;
    double sc_bond_len;
    double sc_steric_radius;
};

struct atom {
    int id;
    bool synthesised;
    char name[MAX_ATOM_NAME_SZ];
    struct vector position, velocity, force;
    double radius;
};

struct residue {
    const struct AA *aa;
    int id;
    bool synthesised;
    size_t num_atoms;
    struct atom *atoms;
};

struct residue *residue_alloc(const struct AA *aa, int id);
void residue_init(struct residue *r, const struct AA *aa, int id);
void residue_free(struct residue *r);
extern struct AA * AA_lookup (register const char *str, register unsigned int len);
void atom_init(struct atom *a, int id, const char *name);

void residue_synth(struct residue *r, struct residue *prev, struct residue *prev2);

#endif /* RESIDUE_H_ */
