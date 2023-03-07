#ifndef TSP_PARSER_H
#define TSP_PARSRE_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Point{
    double x;
    double y;
    int index;
} Point;

typedef struct TSP_data{
    int n_dimensions;
    char* name;
    Point** points;
} TSP_data;

TSP_data parse_file(char* filename);
#endif