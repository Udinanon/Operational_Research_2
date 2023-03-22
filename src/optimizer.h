#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "tsp_data.h"

/// @brief apply multilpe times 2OPT on copy of input TSP Solution
/// @param initial Input tsp solution. is copied and not modified
/// @param data TSP problem used to read information
/// @param iterations number of iterations of 2OPT
/// @return poiter to allocated optimized TSP solution
TSP_solution* two_OPT(TSP_solution* initial, const TSP_data* data, int iterations);

void two_OPT_iteration(TSP_solution* initial, const TSP_data* data);  // apply 2OPT

void swap_array(TSP_solution* solution, int indices[2]);

#endif