#ifndef RESIDUE_H_
#define RESIDUE_H_

struct residue {
    struct AA *aa;
};

struct residue *residue_alloc(const struct AA *aa);
void residue_free(struct residue *r);

#endif /* RESIDUE_H_ */
