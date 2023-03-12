#include "tps_data.h"

TSP_data* init_tsp_data(){
    TSP_data* data = (TSP_data* )malloc(sizeof(TSP_data));
    if(data == NULL){
        logger(FATAL, "TSP_data pointer allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }
    data->n_dimensions = -1;
    data->points = NULL;
    return data;
}

void destroy_tsp_data(TSP_data* data){
    free(data->points);
    free(data);
}

void allocate_points(TSP_data* data){
    if (data->n_dimensions == -1){
        logger(FATAL, "n_dimensions not defined at allocation");
        exit(EXIT_FAILURE);
    }
    data->points = (Point*)malloc(sizeof(Point) * data->n_dimensions);
    if (data->points == NULL) {
        logger(FATAL, "TSP_data points allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }
}