#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "springreader.h"
#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"
#include "bond_angle.h"
#include "rama.h"

#include "cJSON/cJSON.h"

//Macro for perror'ing an error and then goto'ing a label
#define goto_perror(label, error) do { perror(error); goto label; } while(0)
//Similar, but using fprintf rather than perror
#define goto_err(label, ...) do { fprintf(stderr, __VA_ARGS__); goto label; } while(0)
//Similar, but return a value instead of goto'ing
#define ret_err(val, ...) do { fprintf(stderr, __VA_ARGS__); return (val); } while(0)
//Similar, but specifically for cjson errors
#define goto_cjson_error(label, error, json) do {\
    fprintf(stderr, "%s\n", error); \
    print_cjson_error(json); \
    goto label; \
} while(0)

static void print_cjson_error(const char *json_str);
static void set_double_if_set(cJSON *root, const char *key, double *dst);
static void set_int_if_set(cJSON *root, const char *key, int *dst);
static void set_bool_if_set(cJSON *root, const char *key, bool *dst);
static void set_string_if_set(cJSON *root, const char *key, char **dst);


/**
 * Parse  a configuration string into a model.
 *
 * \return A populated model or NULL on error.
 */
struct model * springreader_parse_str(const char *str){
    struct model *m = NULL;
    char *copy = strdup(str);
    if(!copy) goto exit;

    m = model_alloc();
    if(!m) goto free_copy;

    //Read the JSON object
    cJSON *root = cJSON_Parse(str);
    if(!root)
        goto_cjson_error("Error parsing JSON root", str, free_copy);
    //Set config values
    set_double_if_set(root, "timestep", &m->timestep);
    set_double_if_set(root, "synth_time", &m->synth_time);
    set_double_if_set(root, "drag_coefficient", &m->drag_coefficient);
    set_double_if_set(root, "max_synth_angle", &m->max_synth_angle);
    set_double_if_set(root, "until", &m->until);
    set_bool_if_set(root, "use_sterics", &m->use_sterics);
    set_bool_if_set(root, "fix", &m->fix);
    set_bool_if_set(root, "threestate", &m->threestate);
    set_bool_if_set(root, "use_water", &m->use_water);
    set_bool_if_set(root, "shield_drag", &m->shield_drag);
    set_bool_if_set(root, "do_synthesis", &m->do_synthesis);

    read_atoms(root);
    read_springs(root);
    read_angles(root);
    read_torsions(root);

free_copy:
    free(copy);
exit:
    return m;
}

int read_atoms(cJSON *root){
    cJSON *atoms = cJSON_GetObjectItem(root, "atoms");
    if(!atoms)
        ret_err(1, "Couldn't find 'atoms' section\n");
    if(!atoms->child)
        ret_err(1, "No atom records found in 'atoms' section\n");

    cJSON *atom;

    //Go through and determine how many residues we have
    int max_residue = 0;
    int atom_idx = 1;
    for(atom = atoms->child; atom; atom = atom->next, atom_idx++){
        cJSON *res = cJSON_GetObjectItem(atom, "residue");
        if(!res)
            goto_err(error, "Error reading residue from atom %d\n", atom_idx);

        cJSON *res_id = cJSON_GetObjectItem(res, id);
        if(!res_id)
            goto_err(error, "No residue id in atom %d\n", atom_idx);

        if(res_id->valueint > max_residue)
            max_residue = res_id->valueint;
    }

    //Malloc ourselves some residues
    m->num_residues = max_residue;
    m->residues = malloc(sizeof(*m->residues), m->num_residues);
    if(!m->residues)
        goto_perror(error, "Error allocating memory for residues");

    //Now we can get the atoms for each residues
    for(atom = atoms->child; atom; atom = atom->next, atom_idx++){
        int 
    }


error:
    return -1;
}

//Set a value in the model if the key is set in the config file.  For example,
//we want to set drag_coefficient if it is set in the config file, but we want
//to leave the defaults if it is not set.
void set_double_if_set(cJSON *root, const char *key, double *dst){
    cJSON *item = cJSON_GetObjectItem(root, key);
    if(item)
        *dst = item->valuedouble;
}
void set_int_if_set(cJSON *root, const char *key, int *dst){
    cJSON *item = cJSON_GetObjectItem(root, key);
    if(item)
        *dst = item->valueint;
}
void set_bool_if_set(cJSON *root, const char *key, bool *dst){
    cJSON *item = cJSON_GetObjectItem(root, key);
    if(item && item->type == cJSON_False)
        *dst = false;
    else if(item && item->type == cJSON_True)
        *dst = true;
}
void set_string_if_set(cJSON *root, const char *key, char **dst){
    cJSON *item = cJSON_GetObjectItem(root, key);
    if(item)
        *dst = strdup(item->valuestring);
}

//Build a useful error message for cJSON.
void print_cjson_error(const char *json_str){
    int strlen = 0;
    int line = 1;
    char *i;
    const char *error = cJSON_GetErrorPtr();

    //Move error pointer back to start of line (or string)
    while(error > json_str && error[-1] != '\n')
        error--;
    //Find the length of this line
    while(error[strlen] != '\0' && error[strlen] != '\n')
        strlen++;
    //Find the line number we are at
    for(i = index(json_str, '\n'); i < error; i = index(i+1, '\n'))
        line++;

    fprintf(stderr, "Error parsing line %d: %.*s\n", line, strlen, error);
}

struct model * springreader_parse_file(const char *file){
    FILE *f = fopen(file, "r");
    if(!f)
        goto_perror(bail, "Error reading input file");

    //Seek to end of file to get file size
    int rv;
    rv = fseek(f, 0, SEEK_END);
    if(rv == -1)
        goto_perror(close_file, "Error seeking to end of input file");

    long file_sz = ftell(f);
    if(file_sz == -1)
        goto_perror(close_file, "Error getting file size (file too large?)");

    rv = fseek(f, 0, SEEK_SET);
    if(rv == -1)
        goto_perror(close_file, "Error seeking back to start of input file.");

    //Malloc buffer and read file
    char *buf = malloc(file_sz + 1);
    if(!buf)
        goto_perror(close_file, "Error allocating memory to store config file.");

    int nread = fread(buf, 1, file_sz, f);
    if(ferror(f) || nread != file_sz)
        goto_perror(free_buffer, "Error reading from config file");
    //Set terminating nul byte
    buf[nread] = '\0';

    //Close file
    fclose(f);

    struct model *model = springreader_parse_str(buf);
    free(buf);
    return model;

free_buffer:
    free(buf);
close_file:
    fclose(f);
bail:
    return NULL;
}

