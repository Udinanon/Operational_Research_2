#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  
#include <concorde.h>

#include <pthread.h>  

#include <cplex.h>

#define VERBOSE				    10		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)
#define CSVOUT					1

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

//for simulated annealing only
#define TEMPERATURE_COEFF 0.99
                
//data structures  

typedef struct {   
	
	//input data
	int nnodes; 	
	double *xcoord;
	double *ycoord;
	double *costs;							// matrix containing all costs from any node i to any node j

	// parameters 
	int model_type; 
	int method;
	double p;								// Probability to choose the best element
	int len_rcl;							// Restricted Candidate List length
	int refinement;
	int meta;
	int kick;
	int population;							//number of max population for genetic algorithm
	int nrand;
	int randomseed;
	int num_threads;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	int integer_costs;
	int plot;
	// for cplex
	int ncols;

	//global data
	int *succ;								// successors
	int *path;								// path
	int *comp;
	int ncomp;
	double t_start;
	double *best_sol;						// best sol. available    

} instance;        



// chrono.c
double second();

// utils.c
void generate_random_points(double **xcoord, double **ycoord, int nrand, int seed, double x_max, double y_max);
void succ_to_path(int *succ, int *path);
void path_to_succ(int *path, int *succ, int len);
void print_succ(const int *succ, int len);
void print_path(const int *path, int len);
double calculate_succ_cost(int *succ, instance *inst);
double calculate_path_cost(int *path, instance *inst);
int imax(int i1, int i2);
double dmin(double d1, double d2);
double dmax(double d1, double d2);
double random01();

// tsp.c
double cost(int i, int j, instance *inst);

// tsp_meta.c
int two_opt(int *succ, double t_limit, instance *inst);

// tsp_generative.c
int extra_mileage_old(instance *inst);
int extra_mileage(int **out_succ, double p, int len_rcl, double t_limit, instance *inst);
int extra_mileage_compute(int *path, int **out_succ, int len_path, double p, int len_rcl, double t_limit, instance *inst);
int greedy(int **succ, double p, int len_rcl, double t_limit, instance *inst);
int generate_random_path(int **out_path, instance *inst);

// tsp_cplex.c
int TSPopt2(instance *inst);
void addOneCompSec(CPXENVptr env, CPXLPptr lp, int *comp, int comp_number, CPXCALLBACKCONTEXTptr context, instance *inst);
void patchingHeuristicUpdate(CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int nnodes, int ncomp, instance *inst);
void calculateComponents(int **succ, int **comp, int *ncomp, const double *xstar, instance *inst);


// io_utils.c
void plot_best_sol(instance *inst);
void plot_best_sol_old(instance *inst);
void plot_points(instance *inst);
void read_input(instance *inst);
void print_error(const char *err);
void parse_command_line(int argc, char** argv, instance *inst);
void delete_instance(instance *inst);

#endif   /* TSP_H_ */
