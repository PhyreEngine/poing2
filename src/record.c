#include "record.h"
#include "model.h"
#include "vector.h"
#include "residue.h"
#include <stdlib.h>

void record_init(struct record *r, struct model *m, size_t max_records){
    r->max_records = max_records;
    r->natoms = m->num_atoms;

    r->nrecords   = malloc(sizeof(*r->nrecords)   * r->natoms);
    r->avg_jitter = malloc(sizeof(*r->avg_jitter) * r->natoms);
    r->prev_vec   = malloc(sizeof(*r->prev_vec)   * r->natoms);
    r->jitter_buf = malloc(sizeof(*r->jitter_buf) * r->natoms * r->max_records);
    r->buf_idx    = malloc(sizeof(*r->buf_idx)    * r->natoms);

    for(size_t i=0; i < r->natoms; i++){
        r->nrecords[i]   = 0;
        r->avg_jitter[i] = 0;
        r->buf_idx[i]    = 0;
        for(size_t j=0; j < r->max_records; j++)
            r->jitter_buf[i*r->max_records + j] = 0;
    }
}

void record_free(struct record *r){
    free(r->nrecords);
    free(r->avg_jitter);
    free(r->prev_vec);
    free(r->jitter_buf);
    free(r->buf_idx);
}

void record_add(struct record *r, struct model *m){
    for(size_t i=0; i < m->num_atoms; i++){
        if(m->atoms[i].fixed)
            continue;

        if(r->nrecords[i] == 0){
            //If we don't have a reference point, just set that.
            vector_copy_to(&r->prev_vec[i], &m->atoms[i].position);
            r->nrecords[i]++;
        }else{
            //Calculate jitter and copy current position to previous vector
            struct vector displ;
            vsub(&displ, &m->atoms[i].position, &r->prev_vec[i]);
            vector_copy_to(&r->prev_vec[i], &m->atoms[i].position);
            double jitter = vmag(&displ);

            //Calculate the sum of the previous jitters from the current avg.
            //Multiply by (r->nrecords[i] - 1) because of the fencepost problem.
            double sum = r->avg_jitter[i] * (r->nrecords[i] - 1);

            //Index of the circular buffer within the jitter_buf block
            size_t buf_idx = (i * r->max_records);
            //Index within the circular buffer of the item to pop
            size_t curr_idx = r->buf_idx[i];
            size_t next_idx = (r->buf_idx[i] + 1) % r->max_records;

            //If we have reached the maximum number of records to keep, pop the
            //old jitter values. This is the just the next value, because we
            //are using a circular buffer.
            if(r->nrecords[i] == r->max_records)
                sum -= r->jitter_buf[buf_idx + next_idx];
            else
                r->nrecords[i]++;

            r->buf_idx[i] = next_idx;
            r->jitter_buf[buf_idx + curr_idx] = jitter;
            sum += jitter;
            r->avg_jitter[i] = sum / (r->nrecords[i] - 1);
        }
    }
}
