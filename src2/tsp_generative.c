#include "tsp.h"



/*
int extra_mileage_old(instance *inst){

    int *unc_nodes = (int *) calloc(inst->nnodes, sizeof(int));
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    inst->best_sol = (double *) calloc(inst->nnodes, sizeof(double));

    //Choose 2 starting points (max cost)
    double max_cost = 0;
    int temp_i, temp_j, temp_k;
    for(int i=0; i<inst->nnodes; i++){
        for(int j=0; j<i; j++){
            if(cost(i,j,inst) > max_cost){
                max_cost = cost(i,j,inst);
                temp_i = i;
                temp_j = j;
            }
        }
    }
    if(VERBOSE >= 100) printf("Max_cost: %f", max_cost);

    //Add chosen points to successors and set other values to -1
    for(int i=0; i<inst->nnodes; i++) succ[i] = -1;
    succ[temp_i] = temp_j;
    succ[temp_j] = temp_i;
    inst->best_sol[temp_i] = cost(temp_i, succ[temp_i], inst);
    inst->best_sol[temp_j] = cost(temp_j, succ[temp_j], inst);

    //Remove chosen points from uncovered nodes
    int len_unc_nodes = inst->nnodes - 2;
    for(int i=0; i<inst->nnodes; i++){
        unc_nodes[i] = i;
    }
    unc_nodes[temp_i] = unc_nodes[inst->nnodes-1];
    unc_nodes[temp_j] = unc_nodes[inst->nnodes-2];

    double delta, temp_delta;
    int end = 0;
    while(!end){
        delta = INFINITY;
        //for each connection in succ (i,succ[i])
        for(int i=0; i<inst->nnodes; i++){
            if(succ[i] >= 0){
                //for each uncovered node (unc_nodes[j])
                for(int j=0; j<len_unc_nodes; j++){
                    //calculates delta
                    temp_delta = cost(i,unc_nodes[j],inst) + cost(unc_nodes[j],succ[i], inst) - cost(i,succ[i], inst);
                    if(temp_delta < delta){
                        delta = temp_delta;
                        temp_i = i;
                        temp_j = j;
                        temp_k = succ[i];
                    }
                }
            }
        }
        succ[temp_i] = unc_nodes[temp_j];
        succ[unc_nodes[temp_j]] = temp_k;
        inst->best_sol[temp_i] = cost(temp_i, succ[temp_i], inst);
        inst->best_sol[unc_nodes[temp_j]] = cost(unc_nodes[temp_j], succ[unc_nodes[temp_j]], inst);

        if(VERBOSE >= 90) printf("Adding node %d to succ with delta cost: %lf\n", unc_nodes[temp_j], delta);

        unc_nodes[temp_j] = unc_nodes[len_unc_nodes-1];

        len_unc_nodes --;
        if(len_unc_nodes==0){
            end = 1;
        }
    }
    
    //if(VERBOSE >= 70) printf("Len_unc_nodes = %d\n---------------------\n", len_unc_nodes);
    inst->succ = succ;
    if(VERBOSE >= 40){
        printf("Best solution found with EXTRA-MILEAGE algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, (inst->succ)[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("EXTRA-MILEAGE Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
    }

    free(unc_nodes);
    free(succ);
    return 0;
}
*/

int extra_mileage(int **out_succ, double p, int len_rcl, double t_limit, instance *inst){

    //Choose 2 starting points (max cost)
    double max_cost = 0;
    int temp_i, temp_j, temp_k;
    for(int i=0; i<inst->nnodes; i++){
        for(int j=0; j<i; j++){
            if(cost(i,j,inst) > max_cost){
                max_cost = cost(i,j,inst);
                temp_i = i;
                temp_j = j;
            }
        }
    }
    int *starting_path = (int *) calloc(2, sizeof(int));
    starting_path[0] = temp_i;
    starting_path[1] = temp_j;

    extra_mileage_compute(starting_path, &(*out_succ), 2, p, len_rcl, t_limit, inst);
    free(starting_path);
}

// Requires an instance with a path computed and it will
// connect all unconnected nodes
// TODO: pensare a come migliorare il ciclo per non creare nuovamente ogni volta la rcl
int extra_mileage_compute(int *path, int **out_succ, int len_path, double p, int len_rcl, double t_limit, instance *inst){
    
    double t_start = second();
    (*out_succ) = (int *) calloc(inst->nnodes, sizeof(int));
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *rcl = (int *) calloc(len_rcl, sizeof(int));                                    // Restricted Candidate List
    int *rcl_prev = (int *) calloc(len_rcl, sizeof(int));                               // I need to store information of the predecessor
    int *unc_nodes_rcl = (int *) calloc(inst->nnodes-1, sizeof(int));                   // Uncovered nodes flag for Restricted Candidate List
    int *unc_nodes = (int *) calloc(inst->nnodes, sizeof(int));
    


    double best_sol_cost = INFINITY;
    do{

        int len_unc_nodes = inst->nnodes - len_path;
        for(int i=0; i<inst->nnodes; i++) succ[i] = -1;
        path_to_succ(path, succ, len_path);

        int counter = 0;
        for(int i=0; i<inst->nnodes; i++){
            if(succ[i] == -1){
                unc_nodes[counter] = i;
                counter ++;
            }
        }
        if(counter != len_unc_nodes){
            print_error(" error in extra_mileage_compute(), path is shorter than the given len_path");
            return 1;
        }
        // Scan points from uncovered nodes
        int temp_i, temp_j, temp_k;
        double delta, temp_delta;
        int end = 0;
        while(len_unc_nodes > 0){
            // Fill restricted candidate list
            for(int i=0; i<len_unc_nodes; i++) unc_nodes_rcl[i] = 0;
            int end_rcl = 0;
            for(int k=0; (k<len_rcl) && (k<len_unc_nodes); k++){
                delta = INFINITY;
                //for each connection in succ (i,succ[i])
                for(int i=0; i<inst->nnodes; i++){
                    if(succ[i] < 0) continue;
                    //if(succ[i] >= 0){
                        //for each uncovered node    len(unc_nodes[])
                    for(int j=0; j<len_unc_nodes; j++){
                        // skip already visited nodes for RCL.
                        if(unc_nodes_rcl[j] != 0) continue;
                        //calculates delta
                        temp_delta = cost(i,unc_nodes[j],inst) + cost(unc_nodes[j], succ[i], inst) - cost(i, succ[i], inst);
                        if(temp_delta < delta){
                            delta = temp_delta;
                            temp_i = i;
                            temp_j = j;
                            temp_k = succ[i];
                        }
                    }
                    //}
                }
                // Add node with lowest cost to the list
                //printf("Adding node with delta cost: %f\n", delta);
                unc_nodes_rcl[temp_j] = 1;
                rcl_prev[k] = temp_i;
                rcl[k] = temp_j;
                end_rcl ++;

            }
            
            //printf("end_rcl: %d\n", end_rcl);
            if(end_rcl > 0){
                int n = 0;
                if(end_rcl >= 2){
                    if(random01() > p){
                        n = round(random01() * (end_rcl-2)) + 1;
                        // printf("[max %3d] RCL idx: %d, with cost: %f\n", end_rcl, n, cost(rcl_prev[n], unc_nodes[rcl[n]], inst) + cost(succ[rcl_prev[n]], unc_nodes[rcl[n]], inst) 
                        // - cost(rcl_prev[n], succ[rcl_prev[n]], inst));
                    }
                }
                
                succ[unc_nodes[rcl[n]]] = succ[rcl_prev[n]];
                succ[rcl_prev[n]] = unc_nodes[rcl[n]];

                if(VERBOSE >= 90) printf("Adding node %d to succ with delta cost: %lf\n", unc_nodes[rcl[n]], delta);

                unc_nodes[rcl[n]] = unc_nodes[len_unc_nodes-1];

                len_unc_nodes --;
            }else{
                print_error(" unexpected error on function extra_mileage_compute!");
            }
        }
        double sol_cost = calculate_succ_cost(succ, inst);
        if(sol_cost < best_sol_cost){
            for(int i=0; i<inst->nnodes; i++) (*out_succ)[i] = succ[i];
        }
        best_sol_cost = calculate_succ_cost((*out_succ), inst);
        if(VERBOSE >= 50) printf("Solution cost: %f; Best solution cost: %f\n", sol_cost, best_sol_cost);
    }while((second()-t_start) < t_limit);
    
    
    if(VERBOSE >= 40){
        printf("Best solution found with EXTRA-MILEAGE algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, (*out_succ)[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("EXTRA-MILEAGE Total cost: %lf\n", calculate_succ_cost((*out_succ), inst));
    }
    if(CSVOUT >= 1){
        printf("CSVOUT_nodes%d_seed%d;%f;\n", inst->nnodes, inst->randomseed, calculate_succ_cost((*out_succ), inst));
    }


    free(rcl);
    free(rcl_prev);
    free(unc_nodes_rcl);
    free(unc_nodes);
    free(succ);
    return 0;
}

/**
 * @brief GRASP greedy
 * 
 * @param inst problem instance
 * @param p with probability p choose best element; with probability (1-p) choose another element in the RCL(without best element) randomly
 * @param len_rcl length of RCL
 * @return int 1 if any error occurred; 0 otherwise.
 */
int greedy(int **succ, double p, int len_rcl, double t_limit, instance *inst){

    double t_start = second();
    int *unc_nodes = (int *) calloc(inst->nnodes-1, sizeof(int));               // List of uncovered nodes
    (*succ) = (int *) calloc(inst->nnodes, sizeof(int));                        // List of successors
    int *best_sol = (int *) calloc(inst->nnodes, sizeof(int));                  // Best solution
    double best_sol_cost = INFINITY;
    int *rcl = (int *) calloc(len_rcl, sizeof(int));                            // Restricted Candidate List
    int *unc_nodes_rcl = (int *) calloc(inst->nnodes-1, sizeof(int));           // Uncovered nodes flag for Restricted Candidate List

    do{
        for(int i=1; i<inst->nnodes; i++) unc_nodes[i-1] = i;
        for(int i=0; i<inst->nnodes; i++) (*succ)[i] = 0;
        int current_node = 0;
        int next_node;
        double next_node_cost;
        int len_unc_nodes = inst->nnodes - 1;
        while(len_unc_nodes>0){
            // Fill restricted candidate list
            for(int i=0; i<len_unc_nodes; i++) unc_nodes_rcl[i] = 0;
            int end_rcl = 0;
            for(int i=0; (i<len_rcl) && (i<len_unc_nodes); i++){
                next_node_cost = INFINITY;
                for(int j=0; j<len_unc_nodes; j++){

                    if(next_node_cost>cost(current_node, unc_nodes[j], inst) && (unc_nodes_rcl[j] == 0)){
                        next_node = j;
                        next_node_cost = cost(current_node, unc_nodes[j], inst);
                    }
                }
                unc_nodes_rcl[next_node] = 1;
                rcl[i] = next_node;
                end_rcl ++;
            }

            // With probability p decide if choose first element on RCL or another random value in RCL (not first element)
            if(end_rcl == 0) return 1;
            if(end_rcl>=2){
                if(random01() <= p){
                    next_node = rcl[0];
                }else{
                    int n = round(random01() * (end_rcl-2)) + 1;
                    next_node = rcl[n];
                }
            }else{
                next_node = rcl[0];
            }
            (*succ)[current_node] = unc_nodes[next_node];
            current_node = unc_nodes[next_node];

            // Remove the next node from the uncovered nodes
            unc_nodes[next_node] = unc_nodes[len_unc_nodes-1];
            unc_nodes[len_unc_nodes] = -1;
            len_unc_nodes --;

            // At the end add the last connection between first and last nodes.
            if(len_unc_nodes == 0){
                (*succ)[current_node] = 0;
            }

            if(VERBOSE >= 90) printf("Adding node %d to succ with cost: %lf\n", current_node, next_node_cost);
            if(VERBOSE >= 90) printf("Len_unc_nodes = %d\n---------------------\n", len_unc_nodes);
        }


        double sol_cost = calculate_succ_cost((*succ), inst);
        if(sol_cost < best_sol_cost){
            for(int i=0; i<inst->nnodes; i++) best_sol[i] = (*succ)[i];
            best_sol_cost = calculate_succ_cost((*succ), inst);
        }
        if(VERBOSE >= 50) printf("Solution cost: %f; Best solution cost: %f\n", sol_cost, best_sol_cost);


    }while(second()-t_start < t_limit);

    for(int i=0; i<inst->nnodes; i++) (*succ)[i] = best_sol[i];

    if(VERBOSE >= 40){
        printf("Best solution found with GREEDY algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, (*succ)[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("GREEDY Total cost: %lf\n", calculate_succ_cost((*succ), inst));
    }
    if(CSVOUT >= 1){
        printf("CSVOUT_nodes%d_seed%d;%f;\n", inst->nnodes, inst->randomseed, calculate_succ_cost((*succ), inst));
    }

    free(unc_nodes);
    free(best_sol);
    free(rcl);
    free(unc_nodes_rcl);
    return 0;
}



int generate_random_path(int **out_path, instance *inst){
    (*out_path) = (int *) calloc(inst->nnodes, sizeof(int));

    // Start the path from node 0 always.
    int *unc_nodes = (int *) calloc(inst->nnodes-1, sizeof(int));
    int len_unc_nodes = inst->nnodes-1;
    for(int i=0; i<len_unc_nodes; i++) unc_nodes[i] = i+1;

    (*out_path)[0] = 0;
    for(int i=1; i<inst->nnodes; i++){
        //chose the random node among uncovered nodes:
        int n_rand = random01() * (len_unc_nodes-1);
        //add the random node to the path
        (*out_path)[i] = unc_nodes[n_rand];
        //remove the random node from unc_nodes
        unc_nodes[n_rand] = unc_nodes[--len_unc_nodes];
    }
    
    free(unc_nodes);
    return 0;
}







