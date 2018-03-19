#ifndef RATTLE_H_

struct model;
void rattle_push(struct model *m);
void rattle_unconstrained_push(struct model *m);
void rattle_move(struct model *m);

#endif /* RATTLE_H_ */
