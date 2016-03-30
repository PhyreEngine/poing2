#include <stdio.h>
#include <stdlib.h>
#include "../src/linear_spring.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2,
        double epsilon, const char *text){

    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(){
    plan(30);
    struct linear_spring *s;
    struct atom a;
    struct atom b;
    atom_init(&a, 1, "CA");
    atom_init(&b, 1, "CA");
    atom_set_atom_description(&a, atom_description_lookup("CA", 2));
    atom_set_atom_description(&b, atom_description_lookup("CA", 2));

    vector_fill(&a.position, -0.5, 0, 0);
    vector_fill(&b.position, +0.5, 0, 0);
    s = linear_spring_alloc(1, 1.0, &a, &b);

    struct vector result;
    struct vector force1, force2;

    //At equilibrium
    vector_fill(&result, 0, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force = 0 on A");
    is_vector(&force2, &result, 1e-10, "Force = 0 on B");

    //Pushing outwards
    s->distance = 2.0;
    vector_fill(&result, -1, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force = (-1, 0, 0) on A");

    vector_fill(&result, +1, 0, 0);
    is_vector(&force2, &result, 1e-10, "Force = (+1, 0, 0) on B");

    //Pulling inwards
    s->distance = 0.5;
    vector_fill(&result, +0.5, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force = (+0.5, 0, 0) on A");

    vector_fill(&result, -0.5, 0, 0);
    is_vector(&force2, &result, 1e-10, "Force = (-0.5, 0, 0) on B");

    //With a different spring constant
    s->constant = 2.0;
    s->distance = 2.0;
    vector_fill(&result, -2, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force = (-2, 0, 0) on A");

    vector_fill(&result, +2, 0, 0);
    is_vector(&force2, &result, 1e-10, "Force = (+2, 0, 0) on B");

    //Using cutoffs
    s->constant = 1.0;
    s->cutoff   = 0.5;
    s->distance = 1.0;
    vector_fill(&a.position, -1, 0, 0);
    vector_fill(&b.position, +1, 0, 0);
    vector_fill(&result, 0, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force = 0 past cutoff");

    s->cutoff = 2.5;
    vector_fill(&result, 1, 0, 0);
    linear_spring_force(&force1, &force2, s);
    is_vector(&force1, &result, 1e-10, "Force applied within cutoff");

    done_testing();
}
