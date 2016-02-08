#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/bond_angle.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2, 
        double epsilon, const char *text){
    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(int argc, char **argv){
    plan(30);
    struct bond_angle_spring *s;

    struct atom a1, a2, a3;

    vector_fill(&a1.position,  0, 1, 0);
    vector_fill(&a2.position,  0, 0, 0);
    vector_fill(&a3.position,  1, 0, 0);

    s = bond_angle_spring_alloc(&a1, &a2, &a3, 45, 1.0);

    //Try the new force method
    struct vector f1, f2, f3;
    vector_zero(&f1);
    vector_zero(&f2);
    vector_zero(&f3);

    bond_angle_force(&f1, &f2, &f3, s);

    done_testing();
}
