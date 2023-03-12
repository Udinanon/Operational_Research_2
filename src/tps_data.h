#ifndef TSP_DATA_H
#define TSP_DATA_H

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
} TSP_data;

TSP_data* init_tsp_data(); // initialize memory segment on stack

void destroy_tsp_data(TSP_data* ); // deallocate entire structure

void allocate_points(TSP_data* ); // allocate space for points
#endif