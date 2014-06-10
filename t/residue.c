#include <stdio.h>
#include <stdlib.h>
#include "../src/springreader.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

const char *pdb =
"[PDB]\n"
"ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n"
"ATOM      2  CA  GLY A   2       0.000   0.000   0.000\n"
"ATOM      3  CA  GLY A   3       0.000   0.000   0.000\n"
"ATOM      4  CA  GLY A   4       0.000   0.000   0.000\n"
"ATOM      5  CA  GLY A   5       0.000   0.000   0.000\n"
;

int main(){
    plan(10);

    struct model *m = springreader_parse_str(pdb);
    residue_synth(&m->residues[0], NULL, NULL, 3.141);
    residue_synth(&m->residues[1], &m->residues[0], NULL, 3.141);
    residue_synth(&m->residues[2], &m->residues[1], &m->residues[0], 3.141);
    residue_synth(&m->residues[3], &m->residues[2], &m->residues[1], 3.141);
    model_pdb(stdout, m, 0);

    done_testing();
}
