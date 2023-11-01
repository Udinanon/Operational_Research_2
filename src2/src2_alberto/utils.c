#include "tsp.h"



void generate_random_points(double **xcoord, double **ycoord, int nrand, int seed, double x_max, double y_max){
	if(nrand == 0) print_error(" number of random points is not specified");
	srand(2447248+abs(seed));
	for ( int i = 0; i < 1000; i++ ) random01();
	
	*xcoord = (double *) calloc(nrand, sizeof(double)); 	 
	*ycoord = (double *) calloc(nrand, sizeof(double)); 

	for(int i=0; i<nrand; i++){
		(*xcoord)[i] = random01() * x_max;
		(*ycoord)[i] = random01() * y_max;
	}
	
	if(VERBOSE >= 40) printf("Generated %d random points in a rectangle of size %f x %f\n", nrand, x_max, y_max);
}


void succ_to_path(int *succ, int *path){
    int next = 0;
    int counter = 0;
    do{
        path[counter] = next;
        counter ++;
        next = succ[next];
    }while(next!=0);
}


// Writes a successor list given the path(assuming path starts from node 0),
// uncovered nodes in the successor list will mantain their value
void path_to_succ(int *path, int *succ, int len){
    for(int i=0; i<len-1; i++){
        succ[path[i]] = path[i+1];
    }
    succ[path[len-1]] = path[0];
}


void print_succ(const int *succ, int len){
    printf("Successors:\n");
    for(int i=0; i<len; i++){
        printf("{%d: %d} ",i, succ[i]);
    }
    printf("\n");
}

void print_path(const int *path, int len){
    printf("Path:\n");
    int error = 0;
    for(int i=0; i<len; i++){
        printf("{%d: %d} ", i, path[i]);
        for(int j=i+1; j<len; j++){
            if(path[i] == path[j]){
                error = 1;
            }
        }
    }
    printf("\n");
    if(error){
        printf("Error\n");
        exit(1);
    }
}



double calculate_succ_cost(int *succ, instance *inst){
    double total_cost = 0;
    for(int i=0; i<inst->nnodes; i++){
        total_cost += cost(i, succ[i], inst);
    }
    return total_cost;
}

double calculate_path_cost(int *path, instance *inst){
    double total_cost = 0;
    for(int i=0; i<(inst->nnodes-1); i++){
        total_cost += cost(path[i], path[i+1], inst);
    }
    total_cost += cost(path[inst->nnodes-1], path[0], inst);
    return total_cost;
}

int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 
double random01() { return ((double) rand() / RAND_MAX); } // return a random value in range 0.0-1.0

