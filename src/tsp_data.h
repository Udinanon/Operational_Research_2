#ifndef TSP_DATA_H
#define TSP_DATA_H

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "utility.h"

typedef struct Point {
  double x;
  double y;
  int index;
} Point;

typedef struct TSP_data {
  int n_dimensions;
  Point* points;
  double* cost_matrix;
} TSP_data;

typedef struct TSP_solution {
  double cost;
  Point* cycle;
} TSP_solution;

TSP_data* init_tsp_data(); // initialize memory segment on stack

void allocate_points(TSP_data* ); // allocate space for points
void create_cost_matrix(TSP_data* );

void destroy_tsp_data(TSP_data* ); // deallocate entire structure

TSP_solution* NN(TSP_data* ); //Nearest Neighbor

TSP_solution* random_NN(TSP_data* , double ); //NN with randomness

void save_solution(TSP_solution* solution, char* savename); // can't havce circular dependency with Utility so it's here

#endif