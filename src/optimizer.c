#include "optimizer.h"


TSP_solution* two_OPT(TSP_solution* initial, const TSP_data* data, int iterations){
    TSP_solution* solution = allocate_solution(data->n_dimensions + 1);
    solution->cost = initial->cost;
    solution->size = initial->size;
    memcpy(solution->cycle, initial->cycle, sizeof(initial->cycle[0]) * (data->n_dimensions+1));
    for (int i=0; i<iterations; i++){
        two_OPT_iteration(solution, data);
        save_solution(solution, data, "a280", "tmp_sol");
    }
    return solution;
}


void two_OPT_iteration(TSP_solution* solution, const TSP_data* data) {
    int n = solution->size;
    double* cost_matrix = data->cost_matrix;
    double gain = 0;
    int index_to_swap[2] = {-1, -1};
    for (int i = 0; i < n; i++) { // first node of first edge
        for (int j = i + 2; j < n; j++) { // first node of second edge
          double old_cost = cost_matrix[solution->cycle[i + 1] + (solution->cycle[i] * n)] + cost_matrix[solution->cycle[j + 1] + (solution->cycle[j] * n)];
          double new_cost = cost_matrix[solution->cycle[i] + (solution->cycle[j] * n)] + cost_matrix[solution->cycle[j + 1] + (solution->cycle[i + 1] * n)];
          double diff = old_cost - new_cost; // check how muhc we gain
          // gain of swap
          if (diff > gain) {// if advantegeous
            index_to_swap[0] = i;
            index_to_swap[1] = j;
            gain = diff;
            logger(DEBUG, "2OPT found better swap: %d, %d, %f", index_to_swap[0], index_to_swap[1], gain);
          }
        }
  }
  if (gain > 0) {
    logger(INFO, "2OPT best swap: %d, %d, %f", index_to_swap[0], index_to_swap[1], gain);
    swap_array(solution, index_to_swap);
    solution->cost = solution->cost - gain;
  }
}

void swap_array(TSP_solution* solution, int indices[2]){
  //indices[0]+1 is used as we need to swap from i+1 to j including
  int offset = indices[1] - (indices[0]+1);
  int buffer[offset];
  for(int i=0; i<=offset; i++){
    buffer[i] = solution->cycle[indices[1] - i];
  }
  memcpy(&solution->cycle[indices[0] + 1], buffer, offset * sizeof(int));
}
/*
void swap_array(TSP_solution* solution, int indices[2]) {
  int start_index = indices[0] + 1; // we need to swap from i+1 to j to swap the edges (i, i+1) and (j, j+1)
  int end_index = indices[1];
  // (start+end)/2 = half the indeces
  for (int i = 0; i <= ((end_index - start_index) / 2); i++) {
    int temp = solution->cycle[i + start_index];
    solution->cycle[i + start_index] = solution->cycle[end_index - i];
    solution->cycle[end_index - i] = temp;
  }
}*/