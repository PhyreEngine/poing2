#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/torsion_spring.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(vector v1, vector v2, const char *text){
    ok(abs(v1[0] - v2[0]) < 1e-9, "%s: x element", text);
    ok(abs(v1[1] - v2[1]) < 1e-9, "%s: y element", text);
    ok(abs(v1[2] - v2[2]) < 1e-9, "%s: z element", text);
}

int main(int argc, char **argv){
    plan(29);
    struct residue *r1, *r2, *r3, *r4;
    struct torsion_spring *s;

    r1 = residue_alloc(AA_lookup("G", 1));
    r2 = residue_alloc(AA_lookup("G", 1));
    r3 = residue_alloc(AA_lookup("G", 1));
    r4 = residue_alloc(AA_lookup("G", 1));

    r1->position = vector_fill(-1, 0, 0);
    r2->position = vector_fill(0, 0, 0);
    r3->position = vector_fill(0, 0, 1);

    s = torsion_spring_alloc(r1, r2, r3, r4, 90, 1.0);

    //Angle tests - do a few to make sure the geometry is right
    r4->position = vector_fill(1, 0, 1);
    ok(abs(torsion_spring_angle(s) - 180) < 1e-9, "Angle of 180 degrees");

    r4->position = vector_fill(0, 1, 1);
    ok(abs(torsion_spring_angle(s) + 90) < 1e-9, "Angle of -90 degrees");

    r4->position = vector_fill(0, -1, 1);
    ok(abs(torsion_spring_angle(s) - 90) < 1e-9, "Angle of +90 degrees");

    r4->position = vector_fill(1, -1, 1);
    ok(abs(torsion_spring_angle(s) - 135) < 1e-9, "Angle of +135 degrees");

    r4->position = vector_fill(1, 1, 1);
    ok(abs(torsion_spring_angle(s) + 135) < 1e-9, "Angle of -135 degrees");

    //Axis test - should be from r2 to r3
    is_vector(torsion_spring_axis(s), vector_fill(0, 0, 1), "Axis");

    //Get the torque - try with a couple of different angles and constants
    r4->position = vector_fill(-sqrt(2), -sqrt(2), 1);
    s->angle = 0;
    s->constant = 1.0;
    is_vector(
            torsion_spring_torque(s),
            vector_fill(0, 0, -45.0 / 180 * M_PI),
            "Torque at +45 degrees");

    r4->position = vector_fill(-sqrt(2), sqrt(2), 1);
    s->angle = 45;
    s->constant = 2.0;
    is_vector(
            torsion_spring_torque(s),
            vector_fill(0, 0, 2*90. / 180 * M_PI),
            "Torque at -45 degrees (angle=45)");

    //Calculate the force applied to the atom to supply this torque
    r4->position = vector_fill(0, -1, 1);
    s->angle = 0;
    s->constant = 1;
    is_vector(
            torsion_spring_force(s),
            vector_fill(-M_PI/2, 0, 0),
            "Force test 1");

    r4->position = vector_fill(0, -2, 1);
    is_vector(
            torsion_spring_force(s),
            vector_fill(-M_PI/2/2, 0, 0),
            "Force test 2");

    r4->position = vector_fill(0, 2, 1);
    is_vector(
            torsion_spring_force(s),
            vector_fill(-M_PI/2/4, 0, 0),
            "Force test 3");

    r4->position = vector_fill(2, 2, 1);
    is_vector(
            torsion_spring_force(s),
            vector_fill(-M_PI/2/sqrt(8)/2, +M_PI/2/sqrt(8)/2, 0),
            "Force test 4");

    r4->position = vector_fill(2, -2, 1);
    is_vector(
            torsion_spring_force(s),
            vector_fill(-M_PI/2/sqrt(8)/2, -M_PI/2/sqrt(8)/2, 0),
            "Force test 5");

    done_testing();
}
