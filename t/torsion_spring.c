#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/torsion_spring.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2, 
        double epsilon, const char *text){
    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(int argc, char **argv){
    plan(27);
    struct residue *r1, *r2, *r3, *r4;
    struct torsion_spring *s;

    r1 = residue_alloc(AA_lookup("G", 1));
    r2 = residue_alloc(AA_lookup("G", 1));
    r3 = residue_alloc(AA_lookup("G", 1));
    r4 = residue_alloc(AA_lookup("G", 1));

    vector_fill(&r1->position, -1, 0, 0);
    vector_fill(&r2->position,  0, 0, 0);
    vector_fill(&r3->position,  0, 0, 1);

    s = torsion_spring_alloc(r1, r2, r3, r4, 90, 1.0);

    //Angle tests - do a few to make sure the geometry is right
    vector_fill(&r4->position, 1, 0, 1);
    fis(torsion_spring_angle(s), 180, 1e-10, "Angle of 180 degrees");
    ok(abs(torsion_spring_angle(s) - 180) < 1e-9, "Angle of 180 degrees");

    vector_fill(&r4->position, 0, 1, 1);
    fis(torsion_spring_angle(s), -90, 1e-10, "Angle of -90 degrees");

    vector_fill(&r4->position, 0, -1, 1);
    fis(torsion_spring_angle(s), +90, 1e-10, "Angle of +90 degrees");

    vector_fill(&r4->position, 1, -1, 1);
    fis(torsion_spring_angle(s), +135, 1e-10, "Angle of +135 degrees");

    vector_fill(&r4->position, 1, 1, 1);
    fis(torsion_spring_angle(s), -135, 1e-10, "Angle of -135 degrees");

    //Axis test - should be from r2 to r3
    struct vector axis, result, torque, force;
    torsion_spring_axis(&axis, s);
    vector_fill(&result, 0, 0, 1);
    is_vector(&axis, &result, 1e-10, "Axis");

    //Get the torque - try with a couple of different angles and constants
    s->angle = 0;
    s->constant = 1.0;
    vector_fill(&r4->position, -sqrt(2), -sqrt(2), 1);
    vector_fill(&result, 0, 0, -45.0 / 180 * M_PI);
    torsion_spring_torque(&torque, s);
    is_vector(&torque, &result, 1e-10, "Torque at +45 degrees");

    s->angle = 45;
    s->constant = 2.0;
    vector_fill(&r4->position, -sqrt(2), sqrt(2), 1);
    vector_fill(&result, 0, 0, 2*90. / 180 * M_PI);
    torsion_spring_torque(&torque, s);
    is_vector(&torque, &result, 1e-10, "Torque at -45 degrees (angle=45)");

    //Calculate the force applied to the atom to supply this torque
    s->angle = 0;
    s->constant = 1;
    vector_fill(&r4->position, 2, 2, 1);
    vector_fill(&result, -3*M_PI/16, 3*M_PI/16, 0);
    torsion_spring_force(&force, s, R4);
    is_vector(&force, &result, 1e-10, "Force test 1 (on R4)");

    vector_fill(&result, 0, 3*M_PI/4, 0);
    torsion_spring_force(&force, s, R1);
    is_vector(&force, &result, 1e-10, "Force test 1 (on R1)");

    //Move R4 around to check signs
    vector_fill(&r4->position, 2, -2, 1);
    vector_fill(&result, -3*M_PI/16, -3*M_PI/16, 0);
    torsion_spring_force(&force, s, R4);
    is_vector(&force, &result, 1e-10, "Force test 1 (on R4)");

    vector_fill(&result, 0, -3*M_PI/4, 0);
    torsion_spring_force(&force, s, R1);
    is_vector(&force, &result, 1e-10, "Force test 1 (on R1)");

    done_testing();
}
