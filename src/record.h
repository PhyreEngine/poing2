#ifndef RECORD_H_
#define RECORD_H_

#include <stddef.h>
#include <stdbool.h>

struct model;
struct record {
    ///Number of records over which to calculate the moving average.
    size_t max_records;
    ///Number of atoms for which the average is calculated.
    size_t natoms;

    ///Array of length natoms containing the number of jitters that have been
    //recorded for each atom.
    size_t *nrecords;
    ///Array of length natoms containing the average jitter of each atom.
    double *avg_jitter;

    ///Array of previous vectors, used to calculate the displacement.
    struct vector *prev_vec;
    ///Array of old jitters. This is natoms group of circular buffers, each of
    //length max_records.
    double *jitter_buf;
    ///Index in each buffer.
    size_t *buf_idx;
};

void record_init(struct record *r, struct model *m, size_t max_records);
void record_free(struct record *r);
void record_add(struct record *r, struct model *m);

#endif
