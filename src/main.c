#include "cmd_parser.h"
#include "tsp_data.h"
#include "tsp_parser.h"
#include "utility.h"

int main(int argc, char** argv) {
    logger(INFO, "CMD Args counter %d", argc);
    Parameters params = parse(argc, argv);
    srand(params.seed);
    set_logger_level(params.verbosity);
    logger(INFO, "PARAMS %p FILE %s VERBOSE %d THREADS %d", &params, params.filename, params.verbosity, params.n_threads);
    TSP_data* data = init_tsp_data(); 
    parse_file(params.filename, data);
    logger(INFO, "TSP_DATA: Dims %d", data->n_dimensions);
    
    logger(DEBUG, "Point %d: %d %f %f ", data->n_dimensions / 2, data->points[data->n_dimensions / 2].index, data->points[data->n_dimensions / 2].x, data->points[data->n_dimensions / 2].y);
    logger(DEBUG, "Point %d: %d %f %f ", data->n_dimensions - 2, data->points[data->n_dimensions - 2].index, data->points[data->n_dimensions - 2].x, data->points[data->n_dimensions - 2].y);
    create_cost_matrix(data);
    logger(DEBUG, "Distance between the two nodes: %f", data->cost_matrix[data->n_dimensions / 2 + ((data->n_dimensions - 2) * data->n_dimensions)]);

    TSP_solution* solution = NN(data);
    save_solution(solution, "solution");
    TSP_solution* solution = random_NN(data, 0.1);
    save_solution(solution, "solution_rand");

    TSP_solution* sol = random_NN(data, 0.1);

    return 0;
}