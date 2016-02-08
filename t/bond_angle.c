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

int main(){
    plan(9);
    struct bond_angle_spring *s;

    struct atom a1, a2, a3;

    /*
     * a1
     * |
     * |
     * a2----a3
     */
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
    cmp_ok(f1.c[0], ">=", 0, "a1 moving right");
    cmp_ok(f1.c[1], "<=", 0, "a1 moving down");
    fis(f1.c[2], 0, 1e-3, "a1 approximately stationary in z");

    cmp_ok(f2.c[0], "<=", 0, "a2 moving left");
    cmp_ok(f2.c[1], "<=", 0, "a2 moving down");
    fis(f2.c[2], 0, 1e-3, "a2 approximately stationary in z");

    cmp_ok(f3.c[0], "<=", 0, "a3 moving left");
    cmp_ok(f3.c[1], ">=", 0, "a3 moving up");
    fis(f3.c[2], 0, 1e-3, "a3 approximately stationary in z");

    done_testing();
}
