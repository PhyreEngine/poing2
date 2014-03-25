#ifndef SPRINGREADER_H_
#define SPRINGREADER_H_

#include "model.h"

struct model * springreader_parse_str(const char *str);
struct model * springreader_parse_file(const char *file);

#endif /* SPRINGREADER_H_ */
