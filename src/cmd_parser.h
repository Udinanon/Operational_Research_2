#ifndef CMD_PARSER_H
#define CMD_PARSER_H
#include <argp.h>
#include "utility.h"

typedef struct Parameters {
  char* filename;
  enum LOG_LEVEL verbosity; // need to define levels, will have to read on this, now 0-4
  int seed;
  int n_threads;
} Parameters;

Parameters parse(int , char**);

#endif