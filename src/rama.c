#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <ctype.h>
#include "rama.h"
#include "residue.h"

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

/** Read a file containing several grid points. At each grid point, we have the
 * coordinates of the boundary of the nearest Ramachandran region.
 *
 * For the moment, we will just assume 1 degree bins.
 *
 * Upon error, this function returns a non-zero value.
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
        (*points)[phi * 360 + psi].phi = to_phi;
        (*points)[phi * 360 + psi].psi = to_psi;
    }
read_error:
    free(line);
open_error:
    fclose(fin);
    return retval;
}

/**
 * Find a random favoured point. This is probably pretty biased.
 */
void rama_random_init(struct rama_constraint *rama){
    float phi = ((float)random()) / RAND_MAX * 360;
    float psi = ((float)random()) / RAND_MAX * 360;

    //Round to nearest grid point.
    int phi_grid = (int)(phi + 0.5);
    int psi_grid = (int)(psi + 0.5);

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
    if(phi_f < 0)
        phi_f += 180;
    if(psi_f < 0)
        psi_f += 180;

    //Round to nearest grid point.
    int phi_grid = (int)(phi_f + 0.5);
    int psi_grid = (int)(psi_f + 0.5);

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
        //Update the torsion spring.
        rama->phi->angle = closest->phi;
        rama->psi->angle = closest->psi;
        rama->enabled = true;
    }
error:
    return retval;
}

/**
 * Initialise a Ramachandran constraint.
 */
void rama_init(struct rama_constraint *rama,
        struct residue *residue,
        struct residue *next_residue,
        const char *type){

    rama->type = rama_parse_type(type);
    rama->enabled = true;

    struct atom *phi_C  = residue_get_atom(residue, "C");
    struct atom *phi_N  = residue_get_atom(residue, "N");
    struct atom *phi_CA = residue_get_atom(residue, "CA");
    struct atom *phi_next_C = residue_get_atom(next_residue, "C");
    //TODO: Get a good constant
    rama->phi = torsion_spring_alloc(phi_C, phi_N, phi_CA, phi_next_C, 0, 1);

    struct atom *psi_N  = residue_get_atom(residue, "N");
    struct atom *psi_CA = residue_get_atom(residue, "CA");
    struct atom *psi_C  = residue_get_atom(residue, "C");
    struct atom *psi_next_N = residue_get_atom(next_residue, "N");
    rama->psi = torsion_spring_alloc(psi_N, psi_CA, psi_C, psi_next_N, 0, 1);
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
