#include <stdio.h>
#include "cmd_parser.h"
#include "tsp_parser.h"
#include "utility.h"

int main(int argc, char** argv) {
    logger(INFO, "CMD Args counter %d", argc);
    Parameters params = parse(argc, argv);
    set_logger_level(params.verbosity);
    logger(INFO, "PARAMS %p FILE %s VERBOSE %d THREADS %d", &params, params.filename, params.verbosity, params.n_threads);
    TSP_data data = parse_file(params.filename);
    logger(INFO, "TSP_DATA: Dims %d", data.n_dimensions);
    logger(DEBUG, "Point %d: %d %d %d ", data.n_dimensions / 2, data.points[data.n_dimensions / 2]->index, data.points[data.n_dimensions / 2]->x, data.points[data.n_dimensions / 2]->y);
    logger(DEBUG, "Point %d: %d %d %d ", data.n_dimensions - 2, data.points[data.n_dimensions - 2]->index, data.points[data.n_dimensions - 2]->x, data.points[data.n_dimensions - 2]->y);
    return 0;
}