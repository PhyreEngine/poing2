#include <stdio.h>
#include <stdlib.h>
#include "../src/rk4.h"
#include "../src/springreader.h"
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/linear_spring.h"
#include "../src/torsion_spring.h"
#include "../src/vector.h"
#include "tap.h"


int main(int argc, char **argv){
    plan(2);
/*
    struct residue *r1 = residue_alloc(AA_lookup("G", 1));
    struct residue *r2 = residue_alloc(AA_lookup("G", 1));
    struct residue residues[2];
    residues[0] = *r1;
    residues[1] = *r2;

    struct linear_spring *s = linear_spring_alloc(2.0, 1.0, r1, r2);
    vector_fill(r1->position, -0.5, 0, 0);
    vector_fill(r2->position, +0.5, 0, 0);
    struct linear_spring springs[1];
    springs[0] = *s;

    struct model *m = model_alloc();
    m->num_residues        = 2;
    m->num_linear_springs  = 1;
    m->num_torsion_springs = 0;
    m->residues            = residues;
    m->linear_springs      = springs;

    rk4_push(m, 0.1);
    ok(r1->position[0] < -0.5, "Moved r1 in the correct direction");
    ok(r2->position[0] > +0.5, "Moved r2 in the correct direction");

    for(int i=0; i < 1000; i++){
        printf("%g %g %g %g %g %g %g %g %g %g %g %g\n",
                r1->position[0], r1->position[1], r1->position[2],
                r2->position[0], r2->position[1], r2->position[2],
                r1->velocity[0], r1->velocity[1], r1->velocity[2],
                r2->velocity[0], r2->velocity[1], r2->velocity[2]);
        rk4_push(m, 0.1);
    }
    */

    const char *spec = 
        "sequence = GGGG\n"
        "[Linear]\n"
        "1 2 1.0 1.0\n"
        "2 3 1.0 1.0\n"
        "3 4 1.0 1.0\n"
        "[Torsion]\n"
        "1 2 3 4 -45 0.01\n"
        ;

    struct model *m = springreader_parse_str(spec);
    vector_fill(&m->residues[0].position, -1, 0, 0);
    vector_fill(&m->residues[1].position,  0, 0, 0);
    vector_fill(&m->residues[2].position,  0, 0, 1);
    vector_fill(&m->residues[3].position,  0, 1, 1);

    FILE *fout = fopen("test.csv", "w");
    for(int i=0; i < 10000; i++){
        if(i % 10 == 0){
            char *pdb = model_pdb(m);
            fprintf(fout, "MODEL     % d\n", i);
            fprintf(fout, "%s", pdb);
            fprintf(fout, "CONECT% 5d% 5d\n", 1, 2);
            fprintf(fout, "CONECT% 5d% 5d\n", 2, 3);
            fprintf(fout, "CONECT% 5d% 5d\n", 3, 4);
            fprintf(fout, "ENDMDL\n");
            free(pdb);
        }
        rk4_push(m, 0.01);
    }
    fclose(fout);



    done_testing();
}
