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
#define goto_cjson_error(label, json, ...) do {\
    fprintf(stderr, __VA_ARGS__); \
    print_cjson_error(json); \
    goto label; \
} while(0)

static void print_cjson_error(const char *json_str);
static void set_double_if_set(cJSON *root, const char *key, double *dst);
static void set_int_if_set(cJSON *root, const char *key, int *dst);
static void set_bool_if_set(cJSON *root, const char *key, bool *dst);
static void set_string_if_set(cJSON *root, const char *key, char **dst);

static int read_residues(cJSON *root, struct model *m);
static int read_atoms(cJSON *root, struct model *m);
static int read_springs(cJSON *root, struct model *m);
static int read_angles(cJSON *root, struct model *m);
static int read_torsions(cJSON *root, struct model *m);
static int read_rama(cJSON *root, struct model *m);
static int read_constraints(cJSON *root, struct model *m);

static int check_mandatory_keys(cJSON *root, const char **keys, size_t nkeys,
    const char *fmt);

/**
 * Parse  a configuration string into a model.
 *
 * \return A populated model or NULL on error.
 */
struct model * springreader_parse_str(const char *str){
    struct model *m = NULL;
    char *copy = strdup(str);
    if(!copy)
        goto_err(error, "Error allocating copy of json string\n");

    m = model_alloc();
    if(!m)
        goto_err(error, "Error allocating model\n");

    //Read the JSON object
    cJSON *root = cJSON_Parse(str);
    if(!root)
        goto_cjson_error(free_copy, str, "Error parsing JSON root\n");
    //Set config values
    set_double_if_set(root, "timestep", &m->timestep);
    set_double_if_set(root, "synth_time", &m->synth_time);
    set_double_if_set(root, "drag_coefficient", &m->drag_coefficient);
    set_double_if_set(root, "max_synth_angle", &m->max_synth_angle);
    set_double_if_set(root, "until", &m->until);
    set_double_if_set(root, "record_time", &m->record_time);
    set_bool_if_set(root, "use_sterics", &m->use_sterics);
    set_bool_if_set(root, "fix", &m->fix);
    set_bool_if_set(root, "threestate", &m->threestate);
    set_bool_if_set(root, "use_water", &m->use_water);
    set_bool_if_set(root, "shield_drag", &m->shield_drag);
    set_bool_if_set(root, "do_synthesis", &m->do_synthesis);
    set_int_if_set(root, "fix_before", &m->fix_before);

    if(read_residues(root, m))    goto free_copy;
    if(read_atoms(root, m))       goto free_copy;
    if(read_springs(root, m))     goto free_copy;
    if(read_angles(root, m))      goto free_copy;
    if(read_torsions(root, m))    goto free_copy;
    if(read_rama(root, m))        goto free_copy;
    if(read_constraints(root, m)) goto free_copy;

    free(copy);
    return m;

free_copy:
    free(copy);
error:
    free(m);
    return NULL;
}

int read_residues(cJSON *root, struct model *m){
    cJSON *sequence = cJSON_GetObjectItem(root, "sequence");
    if(!sequence)
        ret_err(1, "Couldn't find 'sequence' key\n");
    if(!sequence->valuestring)
        ret_err(1, "The 'sequence' key must be a string\n");

    //Allocate buffer of residues
    int nres = strlen(sequence->valuestring);
    m->residues = malloc(sizeof *(m->residues) * nres);
    m->num_residues = nres;
    if(!m->residues)
        ret_err(1, "Error allocating residues array\n");

    for(size_t i=0; i < nres; i++){
        struct AA *aa = AA_lookup(sequence->valuestring + i, 1);
        if(!aa)
            goto_err(free_res, "Unknown amino acid '%c' at position %lu\n",
                    sequence->valuestring[i], i);

        residue_init(&m->residues[i], i + 1, aa->threeletter);
    }
    return 0;

free_res:
    free(m->residues);
    return 1;
}

int read_atoms(cJSON *root, struct model *m){
    cJSON *atoms = cJSON_GetObjectItem(root, "atoms");
    if(!atoms)
        ret_err(1, "Couldn't find 'atoms' section\n");
    if(!atoms->child)
        ret_err(1, "No atom records found in 'atoms' section\n");

    //Malloc the atoms array
    int natoms = cJSON_GetArraySize(atoms);
    m->num_atoms = natoms;
    m->atoms = malloc(sizeof(*m->atoms) * natoms);
    if(!m->atoms)
        goto_perror(error, "Error allocating atoms array\n");

    //Get all the atoms
    const char *mandatory[3] = {"id", "name", "residue"};
    size_t i = 0;
    for(cJSON *atom = atoms->child; atom; atom = atom->next, i++){
        int fail = check_mandatory_keys(atom, mandatory, 3,
                "Atom no. %d missing key '%s'\n");
        if(fail)
            goto free_atoms;

        cJSON *id      = cJSON_GetObjectItem(atom, "id");
        cJSON *name    = cJSON_GetObjectItem(atom, "name");
        cJSON *residue = cJSON_GetObjectItem(atom, "residue");

        struct atom_description *desc = atom_description_lookup(
                name->valuestring, strlen(name->valuestring));
        if(!desc)
            goto_err(free_atoms, "Unknown atom type '%s'\n", name->valuestring);

        if(id->valueint < 1 || id->valueint > m->num_atoms)
            goto_err(free_atoms, "Atom ID %d out of range at atom %lu\n",
                    id->valueint, i + 1);

        if(residue->valueint < 1 || residue->valueint > m->num_residues)
            goto_err(free_atoms, "Residue index %d out of range at atom %lu\n",
                    residue->valueint, i+1);

        struct atom *a = &m->atoms[id->valueint-1];
        atom_init(a, id->valueint, name->valuestring);
        a->residue_idx = residue->valueint - 1;
        atom_set_atom_description(a, desc);
    }
    return 0;

free_atoms:
    free(m->atoms);
error:
    return -1;
}

int read_constraints(cJSON *root, struct model *m){
    cJSON *constraints = cJSON_GetObjectItem(root, "constraints");

    //It's valid to have no springs, so this isn't an error
    if(!constraints)
        return 0;

    //Find the number of springs we have and malloc an array
    int nsprings = cJSON_GetArraySize(constraints);
    m->num_constraints = nsprings;
    m->constraints = malloc(sizeof(*m->constraints) * nsprings);
    if(!m->constraints)
        goto_perror(alloc_err, "Error allocating constraints array\n");

    const char *mandatory[2] = {"atoms", "distance"};
    size_t i=0;
    for(cJSON *spring = constraints->child; spring; spring = spring->next, i++){
        int fail = check_mandatory_keys(spring, mandatory, 2,
                "Constraint no. %d missing key '%s'\n");
        if(fail)
            goto free_springs;

        cJSON *atoms    = cJSON_GetObjectItem(spring, "atoms");
        cJSON *distance = cJSON_GetObjectItem(spring, "distance");

        if(atoms->type != cJSON_Array)
            goto_err(free_springs,
                    "Key 'atoms' must be an array in spring %lu\n", i+1);
        if(cJSON_GetArraySize(atoms) != 2)
            goto_err(free_springs,
                    "Key 'atoms' must contain 2 atoms in spring %lu\n", i+1);

        for(cJSON *a = atoms->child; a; a = a->next){
            if(a->valueint - 1 < 0 || a->valueint - 1 >= m->num_atoms)
                goto_err(free_springs,
                        "Atom %d does not exist in spring %lu\n",
                        a->valueint, i + 1);
        }

        int a1 = cJSON_GetArrayItem(atoms, 0)->valueint - 1;
        int a2 = cJSON_GetArrayItem(atoms, 1)->valueint - 1;

        struct constraint *s = &m->constraints[i];
        s->a = a1;
        s->b = a2;
        s->distance = distance->valuedouble;
    }
    return 0;

free_springs:
    free(m->constraints);
alloc_err:
    return 1;
}

int read_springs(cJSON *root, struct model *m){
    cJSON *springs = cJSON_GetObjectItem(root, "linear");

    //It's valid to have no springs, so this isn't an error
    if(!springs)
        return 0;

    //Find the number of springs we have and malloc an array
    int nsprings = cJSON_GetArraySize(springs);
    m->num_linear_springs = nsprings;
    m->linear_springs = malloc(sizeof(*m->linear_springs) * nsprings);
    if(!m->linear_springs)
        goto_perror(alloc_err, "Error allocating springs array\n");

    const char *mandatory[2] = {"atoms", "distance"};
    size_t i=0;
    for(cJSON *spring = springs->child; spring; spring = spring->next, i++){
        int fail = check_mandatory_keys(spring, mandatory, 2,
                "Spring no. %d missing key '%s'\n");
        if(fail)
            goto free_springs;

        cJSON *atoms    = cJSON_GetObjectItem(spring, "atoms");
        cJSON *distance = cJSON_GetObjectItem(spring, "distance");
        cJSON *constant = cJSON_GetObjectItem(spring, "constant");
        cJSON *cutoff   = cJSON_GetObjectItem(spring, "cutoff");
        cJSON *hand     = cJSON_GetObjectItem(spring, "handedness");

        if(atoms->type != cJSON_Array)
            goto_err(free_springs,
                    "Key 'atoms' must be an array in spring %lu\n", i+1);
        if(cJSON_GetArraySize(atoms) != 2)
            goto_err(free_springs,
                    "Key 'atoms' must contain 2 atoms in spring %lu\n", i+1);

        for(cJSON *a = atoms->child; a; a = a->next){
            if(a->valueint - 1 < 0 || a->valueint - 1 >= m->num_atoms)
                goto_err(free_springs,
                        "Atom %d does not exist in spring %lu\n",
                        a->valueint, i + 1);
        }

        int a1 = cJSON_GetArrayItem(atoms, 0)->valueint - 1;
        int a2 = cJSON_GetArrayItem(atoms, 1)->valueint - 1;
        double constant_f = (constant)
            ? constant->valuedouble
            : DEFAULT_SPRING_CONSTANT;


        struct linear_spring *s = &m->linear_springs[i];
        linear_spring_init(s, distance->valuedouble, constant_f,
                &m->atoms[a1], &m->atoms[a2]);
        if(cutoff)
            s->cutoff = cutoff->valuedouble;

        if(hand){
            const char *mand_hand[3] = {"inner", "outer", "handedness"};
            fail = check_mandatory_keys(hand, mand_hand, 3,
                    "Spring no. %d missing key '%s'\n");
            if(fail)
                goto free_springs;

            cJSON *inner = cJSON_GetObjectItem(hand, "inner");
            cJSON *outer = cJSON_GetObjectItem(hand, "outer");
            cJSON *handedness = cJSON_GetObjectItem(hand, "handedness");

            //Check atoms exist. Tedious error checking.
            if(inner->valueint - 1 < 0 || inner->valueint - 1 >= m->num_atoms)
                goto_err(free_springs,
                        "Inner atom %d does not exist in spring %lu\n",
                        inner->valueint, i + 1);
            if(outer->valueint - 1 < 0 || outer->valueint - 1 >= m->num_atoms)
                goto_err(free_springs,
                        "Outer atom %d does not exist in spring %lu\n",
                        outer->valueint, i + 1);
            if(strcmp(handedness->valuestring, "RIGHT") != 0
                    && strcmp(handedness->valuestring, "LEFT") != 0)
                goto_err(free_springs,
                        "Handedness of spring %lu is not 'LEFT' or 'RIGHT'\n",
                        i + 1);

            s->inner = &m->atoms[inner->valueint - 1];
            s->outer = &m->atoms[outer->valueint - 1];
            s->right_handed = (strcmp(handedness->valuestring, "RIGHT") == 0);
        }
    }
    return 0;

free_springs:
    free(m->linear_springs);
alloc_err:
    return 1;
}

int read_angles(cJSON *root, struct model *m){
    cJSON *springs = cJSON_GetObjectItem(root, "angle");

    //It's valid to have no springs, so this isn't an error
    if(!springs)
        return 0;

    //Find the number of springs we have and malloc an array
    int nsprings = cJSON_GetArraySize(springs);
    m->num_bond_angles = nsprings;
    m->bond_angles = malloc(sizeof(*m->bond_angles) * nsprings);
    if(!m->bond_angles)
        goto_perror(alloc_err, "Error allocating bond angles array\n");

    const char *mandatory[2] = {"atoms", "angle"};
    size_t i=0;
    for(cJSON *spring = springs->child; spring; spring = spring->next, i++){
        int fail = check_mandatory_keys(spring, mandatory, 2,
                "Angle no. %d missing key '%s'\n");
        if(fail)
            goto free_springs;

        cJSON *atoms    = cJSON_GetObjectItem(spring, "atoms");
        cJSON *angle    = cJSON_GetObjectItem(spring, "angle");
        cJSON *constant = cJSON_GetObjectItem(spring, "constant");
        cJSON *cutoff   = cJSON_GetObjectItem(spring, "cutoff");

        if(atoms->type != cJSON_Array)
            goto_err(free_springs,
                    "Key 'atoms' must be an array in angle %lu\n", i+1);
        if(cJSON_GetArraySize(atoms) != 3)
            goto_err(free_springs,
                    "Key 'atoms' must contain 3 atoms in angle %lu\n", i+1);

        int a1 = cJSON_GetArrayItem(atoms, 0)->valueint - 1;
        int a2 = cJSON_GetArrayItem(atoms, 1)->valueint - 1;
        int a3 = cJSON_GetArrayItem(atoms, 2)->valueint - 1;
        double constant_f = (constant)
            ? constant->valuedouble
            : DEFAULT_BOND_ANGLE_CONST;

        struct bond_angle_spring *s = &m->bond_angles[i];
        bond_angle_spring_init(s, &m->atoms[a1], &m->atoms[a2], &m->atoms[a3],
                angle->valuedouble, constant_f);
        if(cutoff)
            s->cutoff = cutoff->valuedouble;
    }
    return 0;

free_springs:
    free(m->bond_angles);
alloc_err:
    return 1;
}

int read_torsions(cJSON *root, struct model *m){
    cJSON *springs = cJSON_GetObjectItem(root, "torsion");

    //It's valid to have no springs, so this isn't an error
    if(!springs)
        return 0;

    //Find the number of springs we have and malloc an array
    int nsprings = cJSON_GetArraySize(springs);
    m->num_torsion_springs = nsprings;
    m->torsion_springs = malloc(sizeof(*m->torsion_springs) * nsprings);
    if(!m->bond_angles)
        goto_perror(alloc_err, "Error allocating torsion springs array\n");

    const char *mandatory[2] = {"atoms", "angle"};
    size_t i=0;
    for(cJSON *spring = springs->child; spring; spring = spring->next, i++){
        int fail = check_mandatory_keys(spring, mandatory, 2,
                "Torsion constraint no. %d missing key '%s'\n");
        if(fail)
            goto free_springs;

        cJSON *atoms    = cJSON_GetObjectItem(spring, "atoms");
        cJSON *angle    = cJSON_GetObjectItem(spring, "angle");
        cJSON *constant = cJSON_GetObjectItem(spring, "constant");
        cJSON *cutoff   = cJSON_GetObjectItem(spring, "cutoff");

        if(atoms->type != cJSON_Array)
            goto_err(free_springs,
                    "Key 'atoms' must be an array in angle %lu\n", i+1);
        if(cJSON_GetArraySize(atoms) != 4)
            goto_err(free_springs,
                    "Key 'atoms' must contain 4 atoms in angle %lu\n", i+1);

        int a1 = cJSON_GetArrayItem(atoms, 0)->valueint - 1;
        int a2 = cJSON_GetArrayItem(atoms, 1)->valueint - 1;
        int a3 = cJSON_GetArrayItem(atoms, 2)->valueint - 1;
        int a4 = cJSON_GetArrayItem(atoms, 3)->valueint - 1;
        double constant_f = (constant)
            ? constant->valuedouble
            : DEFAULT_BOND_ANGLE_CONST;

        struct torsion_spring *s = &m->torsion_springs[i];
        torsion_spring_init(s,
                &m->atoms[a1], &m->atoms[a2], &m->atoms[a3], &m->atoms[a4],
                angle->valuedouble, constant_f);
        if(cutoff)
            s->cutoff = cutoff->valuedouble;
    }
    return 0;

free_springs:
    free(m->torsion_springs);
alloc_err:
    return 1;
}

int read_rama(cJSON *root, struct model *m){
    cJSON *rama = cJSON_GetObjectItem(root, "ramachandran");

    //It's valid to have no rama constraints, so this isn't an error
    if(!rama)
        return 0;

    const char *top_level[2] = {"data", "constraints"};
    int fail = check_mandatory_keys(rama, top_level, 2,
            "Ramachandran constraint missing key %d %s\n");
    if(fail)
        goto alloc_err;

    //Read file data
    cJSON *data = cJSON_GetObjectItem(rama, "data");
    if(!data)
        goto_err(alloc_err, "No Ramachandran data files given\n");

    for(cJSON *d = data->child; d; d = d->next){
        enum rama_constraint_type type;
        type = rama_parse_type(d->string);
        if(type == UNKNOWN_RAMA)
            goto_err(alloc_err, "Unknown Ramachandran type '%s'\n", d->string);
        if(rama_read_closest(d->valuestring, type))
            goto alloc_err;
    }

    //Find the number of springs we have and malloc an array
    cJSON *springs = cJSON_GetObjectItem(rama, "constraints");
    int nsprings = cJSON_GetArraySize(springs);
    m->num_rama_constraints = nsprings;
    m->rama_constraints = malloc(sizeof(*m->rama_constraints) * nsprings);
    if(!m->rama_constraints)
        goto_perror(alloc_err, "Error allocating Ramachandran array\n");

    const char *mandatory[2] = {"residue", "type"};
    size_t i=0;
    for(cJSON *spring = springs->child; spring; spring = spring->next, i++){
        int fail = check_mandatory_keys(spring, mandatory, 2,
                "Ramachandran cosntraint %d missing key '%s'\n");
        if(fail)
            goto free_springs;

        int residue = cJSON_GetObjectItem(spring, "residue")->valueint;
        const char *type = cJSON_GetObjectItem(spring, "type")->valuestring;

        cJSON *constant = cJSON_GetObjectItem(spring, "constant");
        double constant_f = (constant)
            ? constant->valuedouble
            : DEFAULT_RAMA_CONST;

        if(residue - 1 < 1 || residue >= m->num_residues)
            goto_err(free_springs,
                    "Ramachandran constraint %lu contains out-of-range residues\n",
                    i);

        struct rama_constraint *r = &m->rama_constraints[i];
        rama_init(r, m, residue - 1, type, constant_f);

    }
    return 0;

free_springs:
    free(m->rama_constraints);
alloc_err:
    return 1;
}

int check_mandatory_keys(cJSON *root, const char **keys, size_t nkeys,
    const char *fmt){

    bool failed = false;
    for(size_t i=0; i < nkeys; i++){
        cJSON *value = cJSON_GetObjectItem(root, keys[i]);
        if(!value){
            fprintf(stderr, fmt, i+1, keys[i]);
            failed = true;
        }
    }
    if(failed)
        return 1;
    return 0;
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
    for(i = index(json_str, '\n'); i && i < error; i = index(i+1, '\n'))
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

