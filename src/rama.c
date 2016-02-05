#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <ctype.h>
#include "rama.h"
#include "residue.h"
#include "model.h"

#define NBINS (360*360)

struct phi_psi {
    double phi, psi;
};

struct phi_psi *general     = NULL;
struct phi_psi *alpha       = NULL;
struct phi_psi *beta        = NULL;
struct phi_psi *proline     = NULL;
struct phi_psi *pre_proline = NULL;
struct phi_psi *glycine     = NULL;
struct phi_psi *alanine     = NULL;

static struct phi_psi **rama_data(enum rama_constraint_type type){
    switch(type){
        case GENERAL:     return &general;
        case ALPHA:       return &alpha;
        case BETA:        return &beta;
        case PROLINE:     return &proline;
        case PRE_PROLINE: return &pre_proline;
        case GLYCINE:     return &glycine;
        case ALANINE:     return &alanine;
        default:          return NULL;
    }
}

int rama_is_inited(enum rama_constraint_type type){
    return (*(rama_data(type))) != NULL;
}

/** Read a file containing several grid points. At each grid point, we have the
 * coordinates of the boundary of the nearest Ramachandran region.
 *
 * For the moment, we will just assume 1 degree bins.
 *
 * Upon error, this function returns a non-zero value.
 *
 * This function expects to read values in the range 0--360, not -180--180.
 */
int rama_read_closest(const char *file, enum rama_constraint_type type){
    int retval = 0;

    FILE *fin = fopen(file, "r");
    if(!fin){
        retval = 1;
        error(0, errno, "Error opening %s", file);
        goto open_error;
    }

    struct phi_psi **points = rama_data(type);
    if(!points){
        retval = 1;
        error(0, 0, "Unknown Ramachandran type %d", (int)type);
        goto open_error;
    }
    *points = malloc(sizeof(struct phi_psi) * NBINS);
    if(!(*points)){
        retval = 1;
        error(errno, 0, "Error allocating memory");
        goto open_error;
    }

    //Set angles to -1 for bins inside Ramachandran regions (i.e. not in data)
    for(size_t i=0; i < NBINS; i++){
        (*points)[i].phi = -1;
        (*points)[i].psi = -1;
    }

    int num_points = 0;

    char *line = NULL;
    ssize_t read;
    size_t len = 0;
    while((read = getline(&line, &len, fin)) != -1){
        //Ignore comments
        if(line[0] == '#')
            continue;

        num_points++;

        int phi, psi, to_phi, to_psi;
        int nread = sscanf(line, "%d %d %d %d", &phi, &psi, &to_phi, &to_psi);
        if(nread != 4){
            retval = 1;
            error(0, 0, "Only read %d items from line %d (4 required)",
                    nread, num_points);
            goto read_error;
        }
        //We dont' remove the 180 degree offset when reading so that we can use
        //negative numbers as a flag to indicate missing values.
        (*points)[phi * 360 + psi].phi = to_phi;
        (*points)[phi * 360 + psi].psi = to_psi;
    }
read_error:
    free(line);
    fclose(fin);
open_error:
    return retval;
}

/**
 * Find a random favoured point. This is probably pretty biased.
 */
void rama_random_init(struct rama_constraint *rama){
    float phi = ((float)random()) / RAND_MAX * 360;
    float psi = ((float)random()) / RAND_MAX * 360;

    //Round to nearest grid point.
    int phi_grid = (int)(phi + 0.5) % 360;
    int psi_grid = (int)(psi + 0.5) % 360;

    struct phi_psi **points = rama_data(rama->type);
    struct phi_psi closest = (*points)[phi_grid * 360 + psi_grid];
    if(closest.psi == -1 || closest.psi == -1){
        rama->phi->angle = phi;
        rama->psi->angle = psi;
    }else{
        rama->phi->angle = closest.phi;
        rama->psi->angle = closest.psi;
    }
}

/**
 * Return true if all atoms in this constraint are synthesised.
 */
int rama_is_synthesised(struct rama_constraint *rama){
    return rama->phi->a1->synthesised &&
        rama->phi->a2->synthesised &&
        rama->phi->a3->synthesised &&
        rama->phi->a4->synthesised &&
        rama->psi->a1->synthesised &&
        rama->psi->a2->synthesised &&
        rama->psi->a3->synthesised &&
        rama->psi->a4->synthesised;
}

/**
 * Parse a type string into an enum.
 */
enum rama_constraint_type rama_parse_type(const char *type){
    char copy[strlen(type) + 1];
    for(size_t i=0; i < strlen(type); i++){
        copy[i] = toupper(type[i]);
        copy[i + 1] = '\0';
    }

    if(strcmp(copy, "GENERAL")     == 0) return GENERAL;
    if(strcmp(copy, "ALPHA")       == 0) return ALPHA;
    if(strcmp(copy, "BETA")        == 0) return BETA;
    if(strcmp(copy, "PROLINE")     == 0) return PROLINE;
    if(strcmp(copy, "PRE_PROLINE") == 0) return PRE_PROLINE;
    if(strcmp(copy, "GLYCINE")     == 0) return GLYCINE;
    if(strcmp(copy, "ALANINE")     == 0) return ALANINE;
    return UNKNOWN_RAMA;
};

/**
 * Find the closest Ramachandran region to the given phi/psi angles.
 */
int rama_get_closest(struct rama_constraint *rama){
    int retval = 0;

    double phi_f = torsion_spring_angle(rama->phi);
    double psi_f = torsion_spring_angle(rama->psi);

    //Round to nearest grid point. Remember that the grid goes from 0--360, not
    //-180--180.
    int phi_grid = (int)(phi_f + 180 + 0.5) % 360;
    int psi_grid = (int)(psi_f + 180 + 0.5) % 360;

    struct phi_psi **points = rama_data(rama->type);
    if(!points){
        retval = 1;
        error(0, 0, "Unknown Ramachandran type %d", (int)rama->type);
        goto error;
    }
    if(!(*points)){
        retval = 1;
        error(0, 0, "Ramachandran data not initialised!");
        goto error;
    }

    struct phi_psi *closest = &((*points)[phi_grid * 360 + psi_grid]);


    if(closest->phi < 0 || closest->psi < 0){
        //If we are in an allowed Ramachandran region, disable the constraint.
        rama->enabled = false;
    }else{
        //Update the torsion spring.  Remember that the data files are from
        //0-360, whereas we use -180 -- 180
        rama->phi->angle = closest->phi - 180;
        rama->psi->angle = closest->psi - 180;
        rama->enabled = true;
    }
error:
    return retval;
}

/**
 * Initialise a Ramachandran constraint.
 */
void rama_init(struct rama_constraint *rama,
        const struct model *m,
        size_t residue_idx,
        const char *type,
        float constant){

    rama->type = rama_parse_type(type);
    rama->enabled = true;

    //Find the atoms we want. Doing a linear search is slow, but this only has
    //to be done once for each residue at start up.
    struct atom *phi_prev_C, *phi_N, *phi_CA, *phi_C;
    phi_prev_C = phi_N = phi_CA = phi_C = NULL;

    struct atom *psi_N, *psi_CA, *psi_C, *psi_next_N;
    psi_N = psi_CA = psi_C = psi_next_N = NULL;

    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a = &m->atoms[i];
        if(a->residue_idx == residue_idx - 1){
            if(strcmp(a->name, "C"))
                phi_prev_C = a;
        }else if(a->residue_idx == residue_idx){
            if(strcmp(a->name, "N") == 0){
                phi_N = a;
                psi_N = a;
            }else if(strcmp(a->name, "CA") == 0){
                phi_CA = a;
                psi_CA = a;
            }else if(strcmp(a->name, "C") == 0){
                phi_C = a;
                psi_C = a;
            }
        }else if(a->residue_idx == residue_idx + 1){
            if(strcmp(a->name, "N") == 0)
                psi_next_N = a;
        }
    }
    rama->phi = torsion_spring_alloc(phi_prev_C, phi_N, phi_CA, phi_C, 0, constant);
    rama->psi = torsion_spring_alloc(psi_N, psi_CA, psi_C, psi_next_N, 0, constant);
}

void rama_free_data(){
    free(general);
    free(alpha);
    free(beta);
    free(proline);
    free(pre_proline);
    free(glycine);
    free(alanine);
}
