#ifndef TSP_DATA_H
#define TSP_DATA_H

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

TSP_data* init_tsp_data(); // initialize memory segment on stack

void allocate_points(TSP_data* ); // allocate space for points
void create_cost_matrix(TSP_data* );

void destroy_tsp_data(TSP_data* ); // deallocate entire structure

#endif