#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "residue.h"

struct residue *residue_alloc(const struct AA *aa, int id){
    struct residue *r = malloc(sizeof(struct residue));
    if(!r)
        return NULL;
    residue_init(r, aa, id);
    return r;
}

void residue_init(struct residue *r, const struct AA *aa, int id){
    r->aa = aa;
    r->id = id;
    r->synthesised = false;
    r->num_atoms = 0;
    r->atoms = NULL;
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
    vector_zero(&a->position);
    vector_zero(&a->velocity);
    vector_zero(&a->force);
}

/**
 * Synthesise (i.e. set position and synthesised flag) a new residue.
 *
 * Here, a residue is "synthesised" by placing all atoms into the system and
 * setting the "synthesised" property of the residue and all atoms to true.
 * The residues prev and prev2 are used to define the position of the new atoms
 * in the following way:
 *
 * - If prev and prev2 are NULL, then r is assumed to be the first residue. The
 *   position of the backbone atoms is not changed and the sidechain atoms are
 *   given a random offset with the appropriate bond length.
 *
 * - If only prev2 is NULL, then r is assumed to be the second residue. A
 *   random vector with the mean backbone distance is chosen and the atoms
 *   placed.
 *
 * - If both prev and prev2 are set then a random position is chosen by
 *   starting at the position of the previous backbone atom and choosing a
 *   random vector in the half-sphere about the axis formed by the difference
 *   in the position vectors of the backbone atoms of prev and prev2.
 *
 */
void residue_synth(
        struct residue *r,
        struct residue *prev,
        struct residue *prev2,
        double max_angle){

    struct vector z = {.c = {0, 0, 1}};
    if(!prev2 && ! prev){
        for(size_t i=1; i < r->num_atoms; i++){
            vector_rand(&r->atoms[i].position, M_PI/2-0.1, M_PI/2+0.1);
            vmul_by(&r->atoms[i].position, r->aa->sc_bond_len);
            vadd_to(&r->atoms[i].position, &r->atoms[0].position);
        }
    }else if(prev && !prev2){
        //Put the backbone atom in the Z direction
        vector_rand(&r->atoms[0].position, 0, max_angle / 180 * M_PI);
        vmul_by(&r->atoms[0].position, CA_CA_LEN);
        vadd_to(&r->atoms[0].position, &prev->atoms[0].position);
        //Place any sidechains at a random angle close to the X-Y plane
        for(size_t i=1; i < r->num_atoms; i++){
            vector_rand(&r->atoms[i].position, M_PI/2-0.1, M_PI/2+0.1);
            vmul_by(&r->atoms[i].position, r->aa->sc_bond_len);
            vadd_to(&r->atoms[i].position, &r->atoms[0].position);
        }
    }else if(prev && prev2){
        //Get a random displacment vector for the new CA atom
        vector_rand(&r->atoms[0].position, 0, max_angle / 180 * M_PI);
        vmul_by(&r->atoms[0].position, CA_CA_LEN);
        /*
         *fprintf(stderr, "%g %g %g %g %g %g ", 
         *        prev->atoms[0].position.c[0],
         *        prev->atoms[0].position.c[1],
         *        prev->atoms[0].position.c[2],
         *        r->atoms[0].position.c[0],
         *        r->atoms[0].position.c[1],
         *        r->atoms[0].position.c[2]);
         */

        //Rotate it with the same angle and axis that the (prev - prev2)
        //displacement vector has from the z axis.
        struct vector displ;
        vsub(&displ, &prev->atoms[0].position, &prev2->atoms[0].position);
        /*
         *fprintf(stderr, "%g %g %g ",
         *        displ.c[0],
         *        displ.c[1],
         *        displ.c[2]);
         */

        struct vector rot_axis;
        vcross(&rot_axis, &displ, &z);
        vdiv_by(&rot_axis, vmag(&rot_axis));
        /*
         *fprintf(stderr, "%g %g %g 0 0 1 ",
         *        rot_axis.c[0],
         *        rot_axis.c[1],
         *        rot_axis.c[2]);
         */

        double angle = acos(vdot(&displ, &z) / vmag(&displ));
        vrot_axis(&r->atoms[0].position, &rot_axis, &r->atoms[0].position, -angle);
        /*
         *fprintf(stderr, "%g %g %g\n", 
         *        r->atoms[0].position.c[0],
         *        r->atoms[0].position.c[1],
         *        r->atoms[0].position.c[2]);
         */
        vadd_to(&r->atoms[0].position, &prev->atoms[0].position);

        //Similar procedure with sidechain, but with a random vector starting
        //in the X-Y plane
        for(size_t i=1; i < r->num_atoms; i++){
            vector_rand(&r->atoms[i].position, M_PI/2-0.1, M_PI/2+0.1);
            vmul_by(&r->atoms[i].position, r->aa->sc_bond_len);
            vrot_axis(&r->atoms[i].position, &rot_axis, &r->atoms[i].position, angle);
            vadd_to(&r->atoms[i].position, &r->atoms[0].position);
        }
    }
    r->synthesised = true;
    for(size_t i=0; i < r->num_atoms; i++)
        r->atoms[i].synthesised = true;
}
