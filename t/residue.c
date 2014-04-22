#include <stdio.h>
#include <stdlib.h>
#include "../src/springreader.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

const char *pdb =
"#Sequence:\n"
"sequence = PAPLEQMR\n";

int main(){
    plan(10);

    struct model *m = springreader_parse_str(pdb);
    residue_synth(&m->residues[0], NULL, NULL);
    residue_synth(&m->residues[1], &m->residues[0], NULL);
    residue_synth(&m->residues[2], &m->residues[1], &m->residues[0]);
    residue_synth(&m->residues[3], &m->residues[2], &m->residues[1]);
    model_pdb(stdout, m, 0);

    done_testing();
}
