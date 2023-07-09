#include "tsp.h"

//chrono
double second();
//main
void plot_best_sol(instance *inst);
double random01(); // return a random value in range 0.0-1.0
void print_error(const char *err);

double dist(int i, int j, instance *inst);
void compute_costs(instance *inst);

int comparePaths(int *path1, int *path2, int len);
void copyPaths(int *src, int *dst, int len);


int TSPopt(instance *inst);
int two_opt_tabu(instance *inst, int min_tenure, int max_tenure, double timelimit);
int kick(instance *inst, int n_kick);
int two_opt_vns(instance *inst, int n_kick, double time_limit);
int genetic(instance *inst);
int simulated_annealing(instance *inst);



int TSPopt(instance *inst){
    int error = 0;
    if(inst->model_type == 0){  //greedy search
        if(inst->timelimit == INFINITY){
            if(greedy(&inst->succ, 1, 1, 0, inst)){
                return 1;
            }
        }else{
            if(greedy(&inst->succ, inst->p, inst->len_rcl, inst->timelimit, inst)){
                return 1;
            }
        }

    }
    if(inst->model_type == 1){  //extra mileage
        int path[] = {50, 0};
        //extra_mileage_compute(path, &inst->succ, 2, inst->p, inst->len_rcl, inst->timelimit, inst);
        if(inst->timelimit == INFINITY){
            extra_mileage(&inst->succ, 1, 1, 0, inst);
        }else{
            extra_mileage(&inst->succ, inst->p, inst->len_rcl, inst->timelimit, inst);
        }
    }
    if(inst->model_type == 2) TSPopt2(inst);    //cplex

    if(inst->refinement == 1){
        if(inst->meta == -1){                       //2-OPT algorithm ONLY
            two_opt(inst->succ, INFINITY, inst);
        }else if(inst->meta == 0){                  //TABU search algorithm
            two_opt_tabu(inst, inst->nnodes/10, inst->nnodes/2, inst->timelimit);
        }else if(inst->meta == 1){                  //VNS algorithm
            if(inst->kick != -1){
                two_opt_vns(inst, inst->kick, inst->timelimit);
            }
        }else if(inst->meta == 2){                  //Genetic algorithm
            error = genetic(inst);
            if(error) return error;
        }else if(inst->meta ==3){                   //Simulated annealing
            simulated_annealing(inst);
        }
    }

    return 0;
}



void compute_costs(instance *inst){
    //allocate memory
    inst->costs = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));
    //compute costs
    for(int i=0; i<inst->nnodes; i++){
        inst->costs[i*inst->nnodes + i] = 0;
        for(int j=0; j<i; j++){
            inst->costs[i*inst->nnodes + j] = dist(i,j, inst);
            inst->costs[j*inst->nnodes + i] = inst->costs[i*inst->nnodes + j];
        }
    }
    if(VERBOSE >= 10) printf("Costs computed!\n");
}

double dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.5; // nearest integer 
	return dis+0.0;
}

double cost(int i, int j, instance *inst){
    return inst->costs[i*inst->nnodes + j];
}




int two_opt_tabu(instance *inst, int min_tenure, int max_tenure, double time_limit){
    double t_start = second();
    double best_cost = INFINITY;
    int *best_succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *tabu = (int *) calloc(inst->nnodes, sizeof(int));
    for(int i=0; i<inst->nnodes; i++) tabu[i] = -999999;
    double delta;
    int temp_i, temp_j;
    int tenure;
    
    int iteration = 0;
    int use_max = -1;
    do{
        if(iteration%200 == 0) use_max = -use_max;
        if(use_max == 1){
            tenure = max_tenure;
        }else{
            tenure = min_tenure;
        }
        delta = INFINITY;
        for(int i=0; i<inst->nnodes-1; i++){
            for(int j=i+1; j<inst->nnodes; j++){
                if(inst->succ[i] == j || inst->succ[j] == i) continue;
                if(iteration-tabu[i] <= tenure || iteration-tabu[j] <= tenure) continue;
                double temp_delta = (cost(i,j,inst) + cost(inst->succ[i], inst->succ[j], inst)) - (cost(i, inst->succ[i], inst) + cost(j, inst->succ[j], inst));
                if(temp_delta<delta){
                    delta = temp_delta;
                    temp_i = i;
                    temp_j = j;
                }
            }
        }


        //Change connection between nodes
        int old_succ_i = inst->succ[temp_i];
        int temp = inst->succ[old_succ_i];
        inst->succ[temp_i] = temp_j;
        inst->succ[old_succ_i] = inst->succ[temp_j];

        //Update successors
        int prev = old_succ_i;
        int temp_next = inst->succ[temp];
        while((temp != temp_j) && (temp != temp_next)){
            inst->succ[temp] = prev;
            prev = temp;
            temp = temp_next;
            temp_next = inst->succ[temp];
            if(VERBOSE >=100){
                printf("prev: %d, temp: %d, temp_next: %d\n", prev, temp, temp_next);
                printf("temp_j: %d\n", temp_j);
            }
        }
        inst->succ[temp] = prev;

        if(best_cost > calculate_succ_cost(inst->succ, inst)){
            for(int i=0; i<inst->nnodes; i++) best_succ[i] = inst->succ[i];
            best_cost = calculate_succ_cost(inst->succ, inst);
        }
        if(delta>0){    //worsening move
            tabu[temp_i] = iteration;
            tabu[temp_j] = iteration;
        }
        iteration ++;
        if(VERBOSE >= 90){
            printf("temp_i: %d, temp_j: %d\n",temp_i, temp_j);
            printf("iteration: %d, delta: %f, tenure: %d\n",iteration, delta, tenure);
            printf("Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
        }

    }while(second()-t_start<time_limit);

    for(int i=0; i<inst->nnodes; i++) inst->succ[i] = best_succ[i];

    
    if(VERBOSE >= 40){
        printf("Best solution found with EXTRA-MILEAGE algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->succ[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("EXTRA-MILEAGE Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
    }

    free(best_succ);
    free(tabu);
    return 0;
}

int kick(instance *inst, int n_kick){
    if(VERBOSE>=100) printf("Kicking..\n");
    int *path = (int *) calloc(inst->nnodes, sizeof(int));
    int *new_path = (int *) calloc(inst->nnodes, sizeof(int));
    int *kick_nodes = (int *) calloc(n_kick, sizeof(int));
    int *ordered_nodes = (int *) calloc(n_kick, sizeof(int));
    int prev = 0;
    int counter = 0;
    do{
        prev = inst->succ[prev];
        path[counter] = prev;
        counter ++;
    }while(prev != 0);

    int n_selected = 0;
    while(n_selected<n_kick){
        int n = random01() * (inst->nnodes-1);
        int is_good = 1;
        if((n == 0) || (inst->succ[n] == 0)) is_good = 0;
        for(int i=0; i<n_selected; i++){
            //if the random number was already on the list or was the successor of another number in the list I skip it
            if((n == kick_nodes[i]) || (inst->succ[n] == kick_nodes[i]) || (inst->succ[kick_nodes[i]] == n) ) is_good = 0;
        }
        if(is_good == 1){
            if(VERBOSE >= 100) printf("n: %d, succ[n]: %d\n", n, inst->succ[n]);
            kick_nodes[n_selected] = n;
            n_selected ++;
        }
    }
    
    if(VERBOSE >= 100){ 
        for(int i=0; i<inst->nnodes; i++) printf("%d, ", path[i]);
        printf("\n");
        for(int i=0; i<n_selected; i++) printf("%d, ", kick_nodes[i]);
        printf("End\n");
    }

    //order selected nodes wrt path
    int n_ordered = 0;
    for(int i=0; i<inst->nnodes; i++){
        for(int j=0; j<n_selected; j++){
            if(path[i] == kick_nodes[j]){
                ordered_nodes[n_ordered] = kick_nodes[j];
                kick_nodes[j] = -1;
                n_ordered ++;
            }
        }
    }

    if(VERBOSE >= 100){
        printf("Ordered nodes to connect: ");
        for(int i=0; i<n_ordered; i++) printf("%d, ", ordered_nodes[i]);
        printf("|\n");
    }

    
    //n-kick connections

    //copia tutto fino al primo
    counter = 0;
    int node = path[counter];
    new_path[counter] = node;
    while(node != ordered_nodes[0]){
        counter ++;
        node = path[counter];
        new_path[counter] = node;
        if(VERBOSE >= 100) printf("first_node: %d\n",node);
    }
    //copia dal penultimo all'ultimo dopo il primo
    int counter2 = 0;
    while(path[counter2] != ordered_nodes[n_ordered-2]) counter2 ++;
    while(node != ordered_nodes[n_ordered-1]){

        counter ++;
        counter2 ++;
        node = path[counter2];
        new_path[counter] = node;
        if(VERBOSE >= 100) printf("last_node2: %d\n",node);
    }
    //copia dei blocchi in mezzo
    for(int i=n_ordered-2; i>=1; i--){
        int counter3 = 0;
        while(path[counter3] != ordered_nodes[i-1]) counter3 ++;

        while(node != ordered_nodes[i]){

            counter ++;
            counter3 ++;
            node = path[counter3];
            new_path[counter] = node;
            if(VERBOSE >= 100) printf("node_loop_%d: %d\n",i,node);
        }
    }
    //copia dei rimanenti
    while(node != 0){
        counter ++;
        node = path[counter];
        new_path[counter] = node;
        if(VERBOSE >= 100) printf("node_remaining: %d\n", node);
    }

    if(VERBOSE >= 100){
        printf("New path:\n");
        for(int i=0; i<inst->nnodes; i++) printf("%d, ", new_path[i]);
        printf("\n");
    }

    for(int i=0; i<inst->nnodes-1; i++){
        inst->succ[new_path[i]] = new_path[i+1];
    }
    inst->succ[0] = new_path[0];


    if(VERBOSE >= 60){
        printf("Solution after the kick:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->succ[i]);
        printf("Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
        printf("\n");
    }

    free(path);
    free(new_path);
    free(kick_nodes);
    free(ordered_nodes);
    return 0;
}

int two_opt_vns(instance *inst, int n_kick, double time_limit){
    double t_start = second();
    double best_cost = INFINITY;
    int *best_succ = (int *) calloc(inst->nnodes, sizeof(int));
    double delta;
    int temp_i, temp_j;
    
    double t1;
    do{
        delta = INFINITY;
        for(int i=0; i<inst->nnodes-1; i++){
            for(int j=i+1; j<inst->nnodes; j++){
                if(inst->succ[i] == j || inst->succ[j] == i) continue;
                double temp_delta = (cost(i,j,inst) + cost(inst->succ[i], inst->succ[j], inst)) - (cost(i, inst->succ[i], inst) + cost(j, inst->succ[j], inst));
                if(temp_delta<delta){
                    delta = temp_delta;
                    temp_i = i;
                    temp_j = j;
                }
            }
        }

        if(VERBOSE >= 70){
            printf("temp_i: %d, temp_j: %d\n",temp_i, temp_j);
            printf("delta: %f\n",delta);
        }
        if(delta<0){    //improving move
            //Change connection between nodes
            int old_succ_i = inst->succ[temp_i];
            int temp = inst->succ[old_succ_i];
            inst->succ[temp_i] = temp_j;
            inst->succ[old_succ_i] = inst->succ[temp_j];

            //Update successors
            int prev = old_succ_i;
            int temp_next = inst->succ[temp];
            while((temp != temp_j) && (temp != temp_next)){
                inst->succ[temp] = prev;
                prev = temp;
                temp = temp_next;
                temp_next = inst->succ[temp];
                // printf("prev: %d, temp: %d, temp_next: %d\n", prev, temp, temp_next);
                // printf("temp_j: %d\n", temp_j);
            }
            inst->succ[temp] = prev;
        }else if(delta>0){
            kick(inst, n_kick);
        }
        if(best_cost>calculate_succ_cost(inst->succ, inst)){
            for(int i=0; i<inst->nnodes; i++) best_succ[i] = inst->succ[i];
            best_cost = calculate_succ_cost(inst->succ, inst);
        }

        t1 = second();
    }while(t1-t_start<time_limit);

    for(int i=0; i<inst->nnodes; i++) inst->succ[i] = best_succ[i];

    printf("Best solution found with TWO-OPT + VNS algorithm:\n");
    for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->succ[i]);
    printf("\n");
    printf("Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));

    free(best_succ);
    return 0;
}

int simulated_annealing(instance *inst)
{
    if(inst->timelimit == -1){
        print_error(" time limit not set (use 'tsp -h' for help)");
    }
    double t_start = second();
    double best_cost = INFINITY;
    int* startpath;
    generate_random_path(&startpath, inst);
    double temperature = calculate_path_cost(startpath, inst)/10;
    path_to_succ(startpath, inst->succ, inst->nnodes);
    do{
        int i = random01() * (inst->nnodes - 1);
		int j = random01() * (inst->nnodes - 1);
		while (i == j) j = random01() * (inst->nnodes - 1); //select two random nodes
        double delta = (cost(i,j,inst) + cost(inst->succ[i], inst->succ[j], inst)) - (cost(i, inst->succ[i], inst) + cost(j, inst->succ[j], inst));
        if(VERBOSE >= 80){
            printf("delta: %f\n",delta);
        }
        if(delta<0){//improving move: take the first improving move found
            //Change connection between nodes
            int old_succ_i = inst->succ[i];
            int temp = inst->succ[old_succ_i];
            inst->succ[i] = j;
            inst->succ[old_succ_i] = inst->succ[j];

            //Update successors
            int prev = old_succ_i;
            int temp_next = inst->succ[temp];
            while((temp != j) && (temp != temp_next)){
                inst->succ[temp] = prev;
                prev = temp;
                temp = temp_next;
                temp_next = inst->succ[temp];
                // printf("prev: %d, temp: %d, temp_next: %d\n", prev, temp, temp_next);
                // printf("temp_j: %d\n", temp_j);
            }
            inst->succ[temp] = prev;
            break;
        }
        else{//worsening move
            double prob = exp(-delta/temperature);
            temperature = TEMPERATURE_COEFF*temperature;
            if(random01() <= prob){
                //Change connection between nodes
                int old_succ_i = inst->succ[i];
                int temp = inst->succ[old_succ_i];
                inst->succ[i] = j;
                inst->succ[old_succ_i] = inst->succ[j];

                //Update successors
                int prev = old_succ_i;
                int temp_next = inst->succ[temp];
                while((temp != j) && (temp != temp_next)){
                    inst->succ[temp] = prev;
                    prev = temp;
                    temp = temp_next;
                    temp_next = inst->succ[temp];
                    // printf("prev: %d, temp: %d, temp_next: %d\n", prev, temp, temp_next);
                    // printf("temp_j: %d\n", temp_j);
                }
                inst->succ[temp] = prev;
            }
        }
    }while((second()-t_start) < inst->timelimit);
    succ_to_path(inst->succ, inst->path);
    free(startpath);
    return 0;
}

int genetic(instance *inst){
    if(inst->timelimit == -1){
        print_error(" time limit not set (use 'tsp -h' for help)");
    }
    if(inst->population == -1){
        print_error(" population for genetic algorithm not set (use 'tsp -h' for help)");
    }
    double t_start = second();
    int *population = (int *) calloc(inst->population*inst->nnodes, sizeof(int));
    int *unc_population = (int *) calloc(inst->population, sizeof(int));
    int *parent1 = (int *) calloc((int) inst->population/2, sizeof(int));
    int *parent2 = (int *) calloc((int) inst->population/2, sizeof(int));
    int *new_population = (int *) calloc((int) inst->population/2*inst->nnodes, sizeof(int));
    int best_instance;
    int iteration = 0;
    double best_cost = INFINITY;

    //Generate random population:
    for(int k=0; k<inst->population; k++){
        int *path;
        generate_random_path(&path, inst);
        //copy solution to population
        for(int i=0; i<inst->nnodes; i++){
            population[k*inst->nnodes + i] = path[i];
        }
        if(best_cost>calculate_path_cost(path, inst)){
            best_instance = k;
            best_cost = calculate_path_cost(path, inst);
        }
    }
    printf("Gen: %d, Best solution: %d, with cost: %3.2f\n", iteration, best_instance, best_cost);

    //generate new solutions with genetic algorithm until time ends
    do{
        printf("Iteration: %d, generating new population\n", iteration);
        int len_unc_population = inst->population;
        for(int i=0; i<inst->population; i++) unc_population[i] = i;

        //choose pairs:
        for(int i=0; i<((int) inst->population/2); i++){
            int n = random01() * (len_unc_population-1);
            parent1[i] = unc_population[n];
            unc_population[n] = unc_population[--len_unc_population];

            n = random01() * (len_unc_population-1);
            parent2[i] = unc_population[n];
            unc_population[n] = unc_population[--len_unc_population];
        }
        //print_path(parent1, (int) inst->population/2);
        //print_path(parent2, (int) inst->population/2);

        //generate new instances:
        for(int i=0; i<((int) inst->population/2); i++){

            //copy first part to new population:
            for(int j=0; j<((int) inst->nnodes/2); j++){
                new_population[i*inst->nnodes + j] = population[parent1[i]*inst->nnodes + j];
            }

            //copy second part to new population:
            int counter = (int) inst->nnodes/2;
            for(int j=((int) inst->nnodes/2); j<inst->nnodes; j++){
                int is_good = 1;
                for(int k=0; k<((int) inst->nnodes/2); k++){
                    if(population[parent2[i]*inst->nnodes + j] == population[parent1[i]*inst->nnodes + k]){                //if parent2 has some nodes already visited skip
                        is_good = 0;
                    }
                }
                if(is_good == 1){
                    new_population[i*inst->nnodes + counter] = population[parent2[i]*inst->nnodes + j];
                    counter ++;
                }
            }

            //fix instances using extra-mileage to visit unvisited nodes and then 2OPT algorithm to refine the solution.
            inst->path = &new_population[i*inst->nnodes];
            if(VERBOSE >=100){
                printf("iteration: %d, i: %d\n", iteration, i);
                printf("genetic->generation of new population from: %d, %d:\n", parent1[i], parent2[i]);
                print_path(&population[parent2[i]*inst->nnodes], inst->nnodes);
                print_path(inst->path, counter);
                printf("New instance generated: %d\n", i);
            }
            extra_mileage_compute(inst->path, &inst->succ, counter, 1.0, 1, 0, inst);
            two_opt(inst->succ, INFINITY, inst);
            succ_to_path(inst->succ, &new_population[i*inst->nnodes]);

        }

        //remove equal solutions to increase diversification
        //equal solutions on the new population
        int len_new_population = (int) inst->population/2;
        for(int i=0; i<len_new_population; i++){
            for(int j=i+1; j<len_new_population; j++){
                if(comparePaths(&new_population[i*inst->nnodes], &new_population[j*inst->nnodes], inst->nnodes) == 0){
                    printf("Equal solutions\n");
                    copyPaths(&new_population[(len_new_population-1)*inst->nnodes], &new_population[j*inst->nnodes], inst->nnodes);
                    j = j-1;
                    len_new_population --;
                }
            }
        }
        //equal solutions on the new population wrt previous population
        for(int i=0; i<inst->population; i++){
            for(int j=0; j<len_new_population; j++){
                if(comparePaths(&population[i*inst->nnodes], &new_population[j*inst->nnodes], inst->nnodes) == 0){
                    printf("Equal solutions comparing previous population\n");
                    copyPaths(&new_population[(len_new_population-1)*inst->nnodes], &new_population[j*inst->nnodes], inst->nnodes);
                    j = j-1;
                    len_new_population --;
                }
            }
        }
        
        // update new population
        int n_remove = len_new_population; //(int) inst->population/2;
        int c = 0;
        for(int i=0; i<inst->population; i++){
            if(i == best_instance) continue;          //do not update the best instance of the previous iteration
            if(c == len_new_population) continue;     //scanned all instances on the new_population
            if(calculate_path_cost(&population[i*inst->nnodes], inst) > calculate_path_cost(&new_population[c*inst->nnodes], inst)){
                if(random01() > 0.7){
                    for(int j=0; j<inst->nnodes; j++){
                        population[i*inst->nnodes + j] = new_population[c*inst->nnodes + j];
                    }
                    c ++;
                }
            }
        }
        //update best instance
        for(int k=0; k<inst->population; k++){
            
            if(best_cost>calculate_path_cost(&population[k*inst->nnodes], inst)){
                best_instance = k;
                best_cost = calculate_path_cost(&population[k*inst->nnodes], inst);
            }
        }
        iteration ++;
        printf("Gen: %d, Best solution: %d, with cost: %3.2f\n", iteration, best_instance, best_cost);

    }while((second()-t_start) < inst->timelimit);

    path_to_succ(&population[best_instance*inst->nnodes], inst->succ, inst->nnodes);


    free(population);
    free(unc_population);
    free(parent1);
    free(parent2);
    free(new_population);
    return 0;
}


/**
 * @brief compare 2 paths 
 * 
 * @param path1 first path
 * @param path2 second path
 * @param len length of the paths to compare
 * @return 0 if paths are equals; 1 otherwise.
 */
int comparePaths(int *path1, int *path2, int len){
    for(int i=0; i<len; i++){
        if(path1[i] != path2[i]) return 1;
    }
    return 0;
}

void copyPaths(int *src, int *dst, int len){
    for(int i=0; i<len; i++){
        dst[i] = src[i];
    }
}


