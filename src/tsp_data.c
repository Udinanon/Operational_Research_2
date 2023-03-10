#include "tsp_data.h"

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

TSP_solution* NN(TSP_data* data){
    TSP_solution* sol = (TSP_solution *)malloc(sizeof(TSP_solution));   //allocate the solution
    if(sol == NULL){
        logger(FATAL, "TSP_solution pointer allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }

    int n = data->n_dimensions;

    sol->cycle = (Point*)malloc(sizeof(Point) * (n+1));
    if (sol->cycle == NULL) {
        logger(FATAL, "TSP_solution cycle allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }

    int* already_visited = (int*) calloc(data->n_dimensions, sizeof(int)); //all values to 0
    // initial node is always visited
    Point start = data->points[0];
    sol->cycle[0] = start;
    already_visited[0] = 1;
    double cost = 0;
    int index = 0;

    for (int j = 1; j < n; j++){
        double row[n];
        memcpy(row, &data->cost_matrix[index*n], n*sizeof(double)); //locally copied matrix row

        double min_cost = DBL_MAX;
        int min_pos = -1;
        // check of all nodes whcih are available and find min cost one
        for(int i=0; i<n; i++){
            if(already_visited[i] == 0){
                if (row[i] < min_cost){
                    min_cost = row[i];
                    min_pos = i;
                }
            }
        }

        cost += min_cost;
        already_visited[min_pos] = 1;
        sol->cycle[j] = data->points[min_pos];
        index = min_pos;
    }

    sol->cycle[n] = start;
    cost+=data->cost_matrix[index*n];
    sol->cost = cost;
    return sol;
}

TSP_solution* random_NN(TSP_data* data, double prob){
    
    TSP_solution* sol = (TSP_solution *)malloc(sizeof(TSP_solution));   //allocate the solution
    if(sol == NULL){
        logger(FATAL, "TSP_solution pointer allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }

    int n = data->n_dimensions;

    sol->cycle = (Point*)malloc(sizeof(Point) * (n+1));
    if (sol->cycle == NULL) {
        logger(FATAL, "TSP_solution cycle allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }
    int* already_visited = (int*) calloc(data->n_dimensions, sizeof(int)); //all values to 0
    Point start = data->points[0];
    sol->cycle[0] = start;
    already_visited[0] = 1;
    double cost = 0;
    int index = 0;

    for (int j = 1; j < n; j++){
        double p = get_random();
        double row[n];
        memcpy(row, &data->cost_matrix[index * n], n * sizeof(double));  // locally copied matrix row
        if(p >= prob){
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
            sol->cycle[j] = data->points[min_pos];
            index = min_pos;
        }
        else{
            int available_nodes[n];
            int counter = 0;
            // compile list of available nodes
            for(int i=0; i<n; i++){ 
                if (already_visited[i] == 0){
                    available_nodes[counter] = i;
                    counter++;
                }
            }
            // pick one at random from these
            int random_node = available_nodes[(int)(get_random() * (counter-1))];
            cost += row[random_node];
            already_visited[random_node] = 1;
            sol->cycle[j] = data->points[random_node];
            index = random_node;
        }
    }

    sol->cycle[n] = start;
    cost+=data->cost_matrix[index*n];
    sol->cost = cost;
    return sol;
}

void save_solution(TSP_solution* solution, TSP_data* data, char* problem_name, char* savename) {
    char* filename;
    int size = asprintf(&filename, "./results/%s.txt", savename);
    if (size == -1) {
        logger(FATAL, "Filename allocation failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }
    logger(INFO, "Saving solution to problem %s to filename %s", problem_name, filename);
    FILE* f = fopen(filename, "w");
    if (f == NULL) {
        logger(FATAL, "Save solution file write failed!");
        perror("Error printed by perror");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "%s\n", problem_name);
    for(int i=0; i<data->n_dimensions; i++){
        Point sol_i = solution->cycle[i];
        fprintf(f, "%d ", sol_i.index-1);
    }
    fprintf(f, "\n");
    fprintf(f, "%f\n", solution->cost);
    fclose(f);
    free(filename);
}