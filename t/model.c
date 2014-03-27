#include <stdio.h>
#include <stdlib.h>
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"


int main(int argc, char **argv){
    plan(8);

    struct residue residues[4];
    struct model *m = model_alloc();
    m->timestep = 1;
    m->synth_time = 1;
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

    struct model synth_model;
    model_synth(&synth_model, m);
    cmp_ok(synth_model.num_residues, "==", 0, "0 residues at t=0");

    m->time += 1;
    model_synth(&synth_model, m);
    cmp_ok(synth_model.num_residues, "==", 1, "1 residues at t=1");
    ok(synth_model.residues[0].synthesised,  "Residue 1 synthesised");
    ok(!synth_model.residues[1].synthesised, "Residue 2 not synthesised");
    fis(synth_model.time, m->time, 1e-10, "src attributes copied to dst");

    m->time += 1;
    model_synth(&synth_model, m);
    cmp_ok(synth_model.num_residues, "==", 2, "2 residues at t=2");

    m->time += 1;
    model_synth(&synth_model, m);
    cmp_ok(synth_model.num_residues, "==", 3, "3 residues at t=3");

    done_testing();
}
