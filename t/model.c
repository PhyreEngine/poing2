#include <stdio.h>
#include <stdlib.h>
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"


int main(int argc, char **argv){
    plan(1);

    struct residue residues[4];
    struct model *m = model_alloc();
    m->num_residues = 4;
    m->residues = residues;
    for(size_t i=0; i < m->num_residues; i++)
        m->residues[i].aa = AA_lookup("G", 1);
    vector_fill(&m->residues[0].position, -1, 0, 0);
    vector_fill(&m->residues[1].position,  0, 0, 0);
    vector_fill(&m->residues[2].position,  0, 0, 1);
    vector_fill(&m->residues[3].position,  0, 1, 1);

    const char *pdb_str = model_pdb(m);

    const char *correct =
    "ATOM      1  CA  GLY     1      -1.000   0.000   0.000\n"
    "ATOM      2  CA  GLY     2       0.000   0.000   0.000\n"
    "ATOM      3  CA  GLY     3       0.000   0.000   1.000\n"
    "ATOM      4  CA  GLY     4       0.000   1.000   1.000\n";

    is(pdb_str, correct, "PDB correct");

    done_testing();
}
