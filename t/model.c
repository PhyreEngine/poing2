#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"


int main(){
    plan(8);

    size_t natoms = 20;
    struct residue residues[1];
    struct atom atoms[natoms + 1];

    struct model *m = model_alloc();
    m->timestep = 1;
    m->synth_time = 1;
    m->num_atoms = natoms;
    m->atoms = atoms;
    m->num_residues = 1;
    m->residues = residues;


    residue_init(&residues[0], 1, "ALA");

    for(size_t i=0; i < natoms; i++){
        if(i % 2 == 0){
            atom_init(&m->atoms[i], i+1, "CA");
            atom_set_atom_description(
                    &m->atoms[i],
                    atom_description_lookup("CA", 2));
        }else{
            atom_init(&m->atoms[i], i+1, "O");
            atom_set_atom_description(
                    &m->atoms[i],
                    atom_description_lookup("O", 1));
        }
    }

    struct constraint con;
    con.a = 0;
    con.b = natoms;
    con.distance = 5.0;
    m->num_constraints = 1;
    m->constraints = &con;

    atom_init(&m->atoms[natoms], natoms+1, "O");
    atom_set_atom_description(&m->atoms[natoms], atom_description_lookup("O", 1));
    m->num_atoms++;
    natoms++;

    int npdb = 0;
    for(size_t i=0; i < natoms; i++)
        model_synth_atom(m, i, 20);
    model_pdb(stdout, m, false, &npdb);
    return 0;
}

/**

    m->timestep = 1;
    m->synth_time = 1;
    m->num_residues = 4;
    m->residues = residues;
    for(size_t i=0; i < m->num_residues; i++){
        residue_init(&m->residues[i], i+1);
        strcpy(m->residues[i].name, "GLY");
        m->residues[i].synthesised = true;

        m->residues[i].atoms = malloc(sizeof(struct atom));
        m->residues[i].num_atoms = 1;
        atom_init(&m->residues[i].atoms[0], i+1, "CA");
        m->residues[i].atoms[0].synthesised = true;
    }
    vector_fill(&m->residues[0].atoms[0].position, -1, 0, 0);
    vector_fill(&m->residues[1].atoms[0].position,  0, 0, 0);
    vector_fill(&m->residues[2].atoms[0].position,  0, 0, 1);
    vector_fill(&m->residues[3].atoms[0].position,  0, 1, 1);

    char *pdb_str;
    size_t pdb_sz;
    FILE *pdb_buf = open_memstream(&pdb_str, &pdb_sz);
    model_pdb(pdb_buf, m, false);
    fclose(pdb_buf);

    const char *correct =
    "ATOM      1  CA  GLY     1      -1.000   0.000   0.000\n"
    "ATOM      2  CA  GLY     2       0.000   0.000   0.000\n"
    "ATOM      3  CA  GLY     3       0.000   0.000   1.000\n"
    "ATOM      4  CA  GLY     4       0.000   1.000   1.000\n";

    is(pdb_str, correct, "PDB correct");

    for(size_t i=0; i < m->num_residues; i++)
        m->residues[i].synthesised = false;

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
*/
