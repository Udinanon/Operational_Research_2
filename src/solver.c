#include "solver.h"

TSP_solution* NN(TSP_data* data) {
  return random_NN(data, 0.0);
}

TSP_solution* random_NN(TSP_data* data, double prob) {

  int n = data->n_dimensions;
  TSP_solution* sol = allocate_solution(n+1);

  int* already_visited = (int*)calloc(data->n_dimensions, sizeof(int));  // all values to 0
  int start_index = (int)(get_random() * (data->n_dimensions - 1));
  Point start = data->points[start_index];
  sol->cycle[0] = start.index;
  already_visited[start_index] = 1;
  double cost = 0;
  int index = start_index;

  for (int j = 1; j < n; j++) {
    double p = get_random();
    double row[n];
    memcpy(row, &data->cost_matrix[index * n], n * sizeof(double));  // locally copied matrix row
    if (p >= prob) {
      double min_cost = DBL_MAX;
      int min_pos = -1;
      for (int i = 0; i < n; i++) {
        if (already_visited[i] == 0) {
          if (row[i] < min_cost) {
            min_cost = row[i];
            min_pos = i;
          }
        }
      }
      cost += min_cost;
      already_visited[min_pos] = 1;
      sol->cycle[j] = data->points[min_pos].index;
      index = min_pos;
    } else {
      int available_nodes[n];
      int counter = 0;
      // compile list of available nodes
      for (int i = 0; i < n; i++) {
        if (already_visited[i] == 0) {
          available_nodes[counter] = i;
          counter++;
        }
      }
      // pick one at random from these
      int random_node = available_nodes[(int)(get_random() * (counter - 1))];
      cost += row[random_node];
      already_visited[random_node] = 1;
      sol->cycle[j] = data->points[random_node].index;
      index = random_node;
    }
  }

  sol->cycle[n] = start.index;
  cost += data->cost_matrix[index * n + start_index];
  sol->cost = cost;
  free(already_visited);
  return sol;
}

TSP_solution* Extra_Mileage(TSP_data* data) {
  return Extra_Mileage_partial(data, data->n_dimensions + 1);
}

TSP_solution* Extra_Mileage_partial(TSP_data* data, int ind) {

  int n = data->n_dimensions;
  TSP_solution* sol = allocate_solution(n+1);

  // take the greater distance
  int max = -1;
  int A;
  int B;
  for (int i = 0; i < pow(n, 2); i++) {
    if (data->cost_matrix[i] > max) {
      max = data->cost_matrix[i];
      A = i / n;
      B = i % n;
    }
  }

  int index = 3;
  int* already_visited = (int*)calloc(n + 1, sizeof(int));
  already_visited[0] = A;
  already_visited[1] = B;
  already_visited[2] = A;
  int cost = data->cost_matrix[A * n + B] * 2;

  // find the "least adding cost triangle" for "ind" times (if ind=n+1, I compute a cycle)
  while (index != ind) {
    int available_nodes[n - 2];
    int counter = 0;
    for (int i = 0; i < n; i++) {
      int present = 0;
      for (int j = 0; j < index - 1; j++) {
        if (already_visited[j] == i) {
          present = 1;
          break;
        }
      }
      if (present == 0) {
        available_nodes[counter] = i;
        counter++;
      }
    }
    // in available_nodes there are the indexes of non-cycling nodes; O(n^2) time
    double min = __DBL_MAX__;
    int second;
    int middle;
    for (int i = 0; i < index; i++) {
      for (int j = i + 1; j < index; j++) {  // select the edge i-j
        double cost_ij = data->cost_matrix[already_visited[i] * n + already_visited[j]];
        for (int k = 0; k < counter; k++) {  // select a non-cycling node
          double cost_ik = data->cost_matrix[already_visited[i] * n + available_nodes[k]];
          double cost_kj = data->cost_matrix[available_nodes[k] * n + already_visited[j]];
          double delta = cost_ik + cost_kj - cost_ij;
          if (delta < min) {
            min = delta;
            second = j;
            middle = k;
          }
        }
      }
    }
    cost += min;
    for (int i = index; i > second; i--) {
      already_visited[i] = already_visited[i - 1];
    }
    already_visited[second] = available_nodes[middle];
    index++;
  }
  for (int i = 0; i < n + 1; i++) {
    sol->cycle[i] = data->points[already_visited[i]].index;
  }
  sol->cost = cost;
  free(already_visited);
  return sol;
}