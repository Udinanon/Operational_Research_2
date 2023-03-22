#ifndef SOLVER_H
#define SOLVER_H

#include "tsp_data.h"

TSP_solution* NN(TSP_data*);  // Nearest Neighbor

TSP_solution* random_NN(TSP_data*, double);  // NN with randomness

TSP_solution* Extra_Mileage_partial(TSP_data*, int);

TSP_solution* Extra_Mileage(TSP_data*);  // Extra-Mileage

#endif