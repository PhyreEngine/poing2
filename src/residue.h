#ifndef RESIDUE_H_
#define RESIDUE_H_

#include <stddef.h>
#include <stdbool.h>
#include "vector.h"
#include "model.h"

///Maximum length (including nul byte) of atom name.
#define MAX_ATOM_NAME_SZ 5

///Maximum number of atoms per residue
#define MAX_ATOMS_PER_RES 10

///CA-CA bond length
#define CA_CA_LEN 3.8

///Radius of steric sphere for CA atoms
#define CA_STERIC_RADIUS 1.0

///Structure defining an amino acid type. The definitions of this are supplied
//in AA.gperf.
struct AA {
    const char *oneletter;
    const char *threeletter;
    double mass;
    bool has_sidechain;
    double sc_bond_len;
    double sc_steric_radius;
    double hydrophobicity;
};

struct atom_description {
    const char *name;
    double mass;
    double steric_radius;
    double hydrophobicity;
    bool backbone;
};

struct AA_name {
    const char *oneletter;
    const char *threeletter;
};

struct atom {
    ///1-indexed ID
    int id;
    bool synthesised;
    char name[MAX_ATOM_NAME_SZ];
    struct vector position, velocity, force;
    double radius;
    double mass;
    bool fixed;
    double hydrophobicity;
    size_t residue_idx;
    bool backbone;
};

struct residue {
    int id;
    char name[MAX_ATOM_NAME_SZ];
    bool synthesised;

    size_t num_atoms;
    ///Array of atom IDs
    size_t atoms[MAX_ATOMS_PER_RES];
};

struct residue *residue_alloc(int id, const char *name);
void residue_init(struct residue *r, int id, const char *name);
void residue_free(struct residue *r);
int residue_push_atom(struct residue *r, struct atom *atom);

extern struct AA * AA_lookup (
        register const char *str,
        register unsigned int len);
extern struct atom_description * atom_description_lookup (
        register const char *str,
        register unsigned int len);
void atom_init(struct atom *a, int id, const char *name);
void atom_set_atom_description(struct atom *a,
        const struct atom_description *aa);
struct atom *residue_get_atom(const struct model *m, const struct residue *r,
        const char *name);

void residue_synth(
        struct residue *r,
        struct residue *prev,
        struct residue *prev2,
        double max_angle);

#endif /* RESIDUE_H_ */
