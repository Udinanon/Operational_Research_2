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
    data->cost_matrix = NULL;
    return data;
}

void destroy_tsp_data(TSP_data* data){
    free(data->points);
    free(data->cost_matrix);
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

void create_cost_matrix(TSP_data* data){
    data->cost_matrix = malloc(sizeof(double) * data->n_dimensions * data->n_dimensions); // allocate single large array
    if (data->cost_matrix == NULL) {
        logger(FATAL, "TSP_data cost matrix allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }
    int n = data->n_dimensions;
    double dx = 0, dy = 0;
    Point point_i, point_j;
    for (int i = 0; i < data->n_dimensions; i++){
        for (int j = 0; j < data->n_dimensions; j++){
            if (i==j){
                data->cost_matrix[j+i*n] = 0;
                continue;
            }
            point_i = data->points[i];
            point_j = data->points[j];
            dx = point_i.x -point_j.x;
            dy = point_i.y - point_j.y;
            data->cost_matrix[j + i * n] = sqrt(dx*dx + dy*dy);
        }
    }
}