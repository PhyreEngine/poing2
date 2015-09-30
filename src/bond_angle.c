#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "bond_angle.h"
#include "vector.h"

struct bond_angle_spring * bond_angle_spring_alloc(
        struct atom *a1, struct atom *a2, struct atom *a3,
        double angle, double constant){

    struct bond_angle_spring *spring = malloc(sizeof(struct bond_angle_spring));
    if(!spring)
        return NULL;

    bond_angle_spring_init(spring, a1, a2, a3, angle, constant);
    return spring;
}

void bond_angle_spring_init(
        struct bond_angle_spring *s,
        struct atom *a1, struct atom *a2, struct atom *a3,
        double angle, double constant){
    s->a1 = a1;
    s->a2 = a2;
    s->a3 = a3;
    s->angle = angle;
    s->constant = constant;
    s->enabled = true;
}

void bond_angle_spring_free(struct bond_angle_spring *s){
    free(s);
}

void bond_angle_force(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct bond_angle_spring *s){

    //Let's try a harmonic bond potential. See the GROMACS manual, section
    //4.2.5.

    //Bond vectors (j is atom 2, the central atom)
    struct vector r_ij, r_kj;
    vsub(&r_ij, &s->a1->position, &s->a2->position);
    vsub(&r_kj, &s->a3->position, &s->a2->position);

    //Calculate modulus of bond vectors
    float r_ij_mod = vmag(&r_ij);
    float r_kj_mod = vmag(&r_kj);

    //Calculate theta
    float cos_theta = vdot(&r_ij, &r_kj) / (r_ij_mod * r_kj_mod);

    //Fix floating point inaccuracy
    if(cos_theta < -1)
        cos_theta = 1;
    else if(cos_theta > 1)
        cos_theta = 1;

    float theta = acos(cos_theta);

    if(1.0f - cos_theta*cos_theta == 0){
        vector_zero(f1);
        vector_zero(f2);
        vector_zero(f3);
        return;
    }

    //Calculate d/dtheta part:
    float constant = -s->constant * (theta - (s->angle/180*M_PI));
    //Add the d(acos)/d(cos) part:
    constant *= (-1.0f) / sqrt(1.0f - cos_theta*cos_theta);

    //Now the d(cos theta)/d(r_i) part. This is a vector.
    struct vector r_kj_by_rijkj, r_ij_by_rijij;
    vmul(&r_kj_by_rijkj, &r_kj, 1.0f / (r_ij_mod * r_kj_mod));
    vmul(&r_ij_by_rijij, &r_ij, cos_theta / (r_ij_mod * r_ij_mod));

    struct vector dcostheta_dri;
    vsub(&dcostheta_dri, &r_kj_by_rijkj, &r_ij_by_rijij);

    //And d(cos theta)/d(r_k)
    struct vector r_ij_by_rijkj, r_kj_by_rkjkj;
    vmul(&r_ij_by_rijkj, &r_ij, 1.0f / (r_ij_mod * r_kj_mod));
    vmul(&r_kj_by_rkjkj, &r_kj, cos_theta / (r_kj_mod * r_kj_mod));

    struct vector dcostheta_drk;
    vsub(&dcostheta_drk, &r_ij_by_rijkj, &r_kj_by_rkjkj);

    //Combine parts
    vmul(f1, &dcostheta_dri, constant);
    vmul(f3, &dcostheta_drk, constant);
    vadd(f2, f1, f3);
    vmul_by(f2, -1);
}
