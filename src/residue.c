#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "residue.h"
#include "model.h"

struct residue *residue_alloc(int id, const char *name){
    struct residue *r = malloc(sizeof(struct residue));
    if(!r)
        return NULL;
    residue_init(r, id, name);
    return r;
}

void residue_init(struct residue *r, int id, const char *name){
    r->id = id;
    r->synthesised = false;
    r->num_atoms = 0;
    strncpy(r->name, name, MAX_ATOM_NAME_SZ);
}

void residue_free(struct residue *r){
    free(r);
}

void atom_init(struct atom *a, int id, const char *name){
    a->id = id;
    a->synthesised = false;
    a->name[0] = '\0';
    strncat(a->name, name, MAX_ATOM_NAME_SZ-1);
    a->radius = 0;
    a->fixed = false;
    vector_zero(&a->position);
    vector_zero(&a->velocity);
    vector_zero(&a->force);
    a->hydrophobicity = 0.0;
    a->residue_idx = 0;
}

void atom_set_atom_description(struct atom *a,
        const struct atom_description *aa){
    a->radius = aa->steric_radius;
    a->mass   = aa->mass;
    a->hydrophobicity = aa->hydrophobicity;
    a->backbone = aa->backbone;
}

/**
 * Push an atom onto a residue.
 *
 * This reallocates the atoms array.
 *
 * \param r     Residue onto which to add the atom.
 * \param atom  Atom to add onto the residue.
 *
 * \return Non-zero if adding the atom failed.
 */
int residue_push_atom(struct residue *r, struct atom *atom){
    if(r->num_atoms == MAX_ATOMS_PER_RES)
        return 1;

    r->num_atoms++;
    r->atoms[r->num_atoms - 1] = atom->id;
    return 0;
}

/**
 * Get an atom with the given name from the residue. Returns NULL if the atom is
 * not found.
 */
struct atom *residue_get_atom(const struct model *m, const struct residue *r,
        const char *name){

    for(size_t i=0; i < r->num_atoms; i++){
        struct atom *a = &m->atoms[r->atoms[i]];
        if(strcmp(a->name, name) == 0)
            return a;
    }
    return NULL;
}

