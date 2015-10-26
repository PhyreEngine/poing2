#ifndef LEAPFROG_H_
#define LEAPFROG_H_

#include "model.h"

void leapfrog_init(struct model *model);
void leapfrog_push(struct model *model);

#endif /* LEAPFROG_H_ */

