#include <config.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "springreader.h"
#include "model.h"
#include "rk4.h"
#include "sterics.h"

#ifdef _GNU_SOURCE
#   include <fenv.h>
#endif

static struct option opts[] = { {"help",     no_argument,       0, 'h'},
    {"snapshot",   required_argument, 0, 's'},
    {"no-connect", no_argument,       0, 'c'},
    {"seed",       required_argument, 0, 'r'},
    {"kinetic",    required_argument, 0, 'k'},
    {0, 0, 0, 0}
};
const char *opt_str = "hs:u:s:k:r:";

const char *usage_str =
"Usage: poing [OPTIONS] <SPEC>\n"
"\n"
"Argument SPEC is mandatatory and must be a specification file for a "
"simulation.\n"
"Available options:\n"
"  -h, --help         Display this help message.\n"
"  -s, --snapshot=N   Write a PDB snapshot every N steps.\n"
"  -r, --seed=S       Use fixed random seed S.\n"
"  -k, --kinetic=F    Write kinetic energies to file F.\n"
"      --no-connect   Do not print CONECT records for each spring.\n"
;
int snapshot = -1;
bool print_connect = true;
bool fixed_seed = false;
int random_seed = 0;
char *kinetic = NULL;

void usage(const char *msg, int exitval){
    FILE *out = (exitval < 2) ? stdout : stderr;
    if(msg)
        fprintf(out, "%s\n", msg);
    fprintf(out, "%s", usage_str);
    exit(exitval);
}

char * get_options(int argc, char **argv){
    int c;
    int option_index;
    while((c = getopt_long(argc, argv, opt_str, opts, &option_index)) != -1){
        switch(c){
            case 'h':
                usage(NULL, 1);
            case 's':
                snapshot = atoi(optarg);
                break;
            case 'c':
                print_connect = false;
                break;
            case 'r':
                fixed_seed = true;
                random_seed = atoi(optarg);
                break;
            case 'k':
                kinetic = optarg;
                break;
        }
    }
    if(optind >= argc)
        usage("No specification file supplied.", 2);
    return argv[optind];
}

int main(int argc, char **argv){
#if defined(_GNU_SOURCE) && !defined(__FAST_MATH__)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    if(!fixed_seed){
        srand(time(NULL) * getpid());
    }else{
        srand(random_seed);
    }


    /* Get options and whatnot */
    char * spec = get_options(argc, argv);
    struct model *model = springreader_parse_file(spec);
    if(!model)
        return 2;

    /* Open output file for kinetic energy if specified. */
    FILE *kinetic_out;
    if(kinetic){
        kinetic_out = fopen(kinetic, "w");
        if(!kinetic_out){
            perror("Couldn't open kinetic energy output file");
            exit(1);
        }
    }


    /* Initialise steric grid if it is being used. */
    struct steric_grid *steric_grid = NULL;
    if(model->use_sterics || model->use_water || model->shield_drag){
         steric_grid = steric_grid_alloc(6);
        model->steric_grid = steric_grid;
    }


    /* Run. */
    struct model state;
    unsigned long i=0;
    while(model->time < model->until){

        /* Do synth steps */
        model_synth(&state, model);
        leapfrog_push(&state);

        /* If we're doing a three-state synthesis then fix/unfix residues
         * and enable/disable torsion springs. */
        if(model->threestate && state.num_residues > 0){
            int num_syn  = state.time / state.synth_time;
            double dt    = state.timestep;
            double since = state.time - (num_syn * state.synth_time);
            struct residue *res  = &state.residues[state.num_residues - 1];
            struct residue *prev = NULL;
            if(state.num_residues > 1)
                prev = &state.residues[state.num_residues - 2];

            /* Fix all atoms for the first 2/3 synth_time. For the first
             * synth_time / 3 use linear springs and not torsion springs.
             * Then, for the next synth_time / 3 use torsional springs
             * only. Then enable them all and unfix atoms.
             *
             * This should get the correct distances but not necessarily
             * the correct handedness. Then we fix the handedness with the
             * linear springs turned off so that the torsion spring does
             * not have to overcome the linear springs. All atoms except
             * the curent one are fixed so that effects don't carry on to
             * the rest of the atoms.*/

            double off_lin = 1 * state.synth_time / 3;
            double unfix   = 2 * state.synth_time / 3;

            /* Fix all residues. */
            if(since - dt <= dt){
                fprintf(stderr, "Fixing previous residues and disabling torsion springs\n");

                for(size_t j=0; j < state.num_residues - 1; j++)
                    for(size_t k=0; k < state.residues[j].num_atoms; k++){
                        state.residues[j].atoms[k].fixed = true;
                        vector_zero(&state.residues[j].atoms[k].velocity);
                    }

                for(size_t j=0; j < state.num_torsion_springs; j++)
                    state.torsion_springs[j].enabled = false;

            }else if(since - dt < off_lin && since >= off_lin){
                fprintf(stderr, "Enabling torsion springs, disabling linear springs\n");

                for(size_t j=0; j < state.num_torsion_springs; j++)
                    state.torsion_springs[j].enabled = true;

                for(size_t j=0; j < state.num_linear_springs; j++)
                    state.linear_springs[j].enabled = false;

                if(prev){
                    for(size_t j = 0; j < state.num_linear_springs; j++){
                        //Check if the spring links this residue with the
                        //previous residue.
                        struct linear_spring *spr = &state.linear_springs[j];
                        bool use = false;
                        for(size_t k=0; k < res->num_atoms; k++){
                            for(size_t l=0; l < prev->num_atoms; l++){
                                struct atom *a1 = &res->atoms[k];
                                struct atom *a2 = &prev->atoms[l];
                                if(spr->a->id == a1->id && spr->b->id == a2->id){
                                    use = true;
                                    break;
                                }
                            }
                            if(use)
                                break;
                        }
                        if(use)
                            spr->enabled = true;
                    }
                }

            }else if(since - dt < unfix && since >= unfix){
                fprintf(stderr, "Unfixing atoms, enabling linear springs\n");

                for(size_t j=0; j < state.num_residues; j++)
                    for(size_t k=0; k < state.residues[j].num_atoms; k++)
                        state.residues[j].atoms[k].fixed = false;

                for(size_t j=0; j < state.num_linear_springs; j++)
                    state.linear_springs[j].enabled = true;

            }
        }



        /* Print snapshot at the correct times. */
        if(snapshot > 0 &&
                (int)(state.time / snapshot) > (int)(model->time / snapshot)){
            printf("MODEL     %lu\n", ++i);
            model_pdb(stdout, &state, print_connect);
            printf("ENDMDL\n");

            /* Also print KE if the output file was supplied. */
            if(kinetic){
                double ke = 0;
                for(size_t j=0; j < state.num_residues; j++){
                    for(size_t k=0; k < state.residues[j].num_atoms; k++){
                        struct atom *atom = &state.residues[j].atoms[k];
                        ke += 0.5* atom->mass * vdot(
                                &atom->velocity, &atom->velocity);
                    }
                }
                fprintf(kinetic_out, "%f\t%f\n", state.time, ke);
                fflush(kinetic_out);
            }
        }
        model->time = state.time;
    }
    if(kinetic)
        fclose(kinetic_out);
    steric_grid_free(steric_grid);
    model_free(model);
    return 0;
}
