#include <stdio.h>
#include "cmd_parser.h"
#include "tsp_parser.h"
#include "utility.h"

int main(int argc, char** argv) {
    logger(NULL, INFO, "CMD Args %s counter %d", argv, argc);
    Parameters params = parse(argc, argv);
    logger(NULL, INFO, "PARAMS %p FILE %s VERBOSE %d THREADS %d \n", &params, params.filename, params.verbosity, params.n_threads);
    TSP_data data = parse_file(params.filename);
    return 0;
}