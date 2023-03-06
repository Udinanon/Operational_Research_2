#include <stdio.h>
#include "cmd_parser.h"

int main(int argc, char** argv) {
    Parameters params = parse(argc, argv);
    printf("PARAMS %p FILE %s VERBOSE %d THREADS %d \n", &params, params.filename, params.verbosity, params.n_threads);
    return 0;
}