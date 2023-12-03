#include "tsp.h"

static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) ;
double dist(int i, int j, instance *inst);
void print_error(const char *err);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void localBranchingAddRow(CPXENVptr env, CPXLPptr lp, const double *xstar, int k, instance *inst);
void fixLowerBound(CPXENVptr env, CPXLPptr lp, int k, int succ_k, instance *inst);
void freeLowerBound(CPXENVptr env, CPXLPptr lp, int min_k, int max_k, instance *inst);
void freeAllLowerBound(CPXENVptr env, CPXLPptr lp, instance *inst);
void updateModel(instance *inst, CPXENVptr env, CPXLPptr lp, int ncomp, int* comp);
void patchingHeuristicUpdate(CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int nnodes, int *ncomp, instance *inst, CPXCALLBACKCONTEXTptr);
int xpos(int i, int j, instance *inst);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int TSPopt2(instance *inst);	//he expects a singular function with different number of parameters for Benders: one without vns, one with vns, ecc.
//don't try smarter implementation of benders, like last year

/**************************************************************************************************************************/
int TSPopt2(instance *inst)
/**************************************************************************************************************************/
{  

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0); 
	// ...


	// use the optimal solution found by CPLEX
	
	int ncols = CPXgetnumcols(env, lp);
	inst->ncols = ncols;

	double *xstar;


	if(strncmp(inst->method, "Basic", 5) == 0){											// Benders' Loop method
		do{
			error = CPXmipopt(env,lp);
			if ( error ) 
			{
				printf("CPX error code %d\n", error);
				print_error("CPXmipopt() error"); 
			}
			xstar = (double *) calloc(ncols, sizeof(double));
			if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("CPXgetx() error");	
			calculateComponents(&inst->succ, &inst->comp, &inst->ncomp, xstar, inst);
			updateModel(inst, env, lp, inst->ncomp, inst->comp);
		}while(inst->ncomp > 1);
	}else if(strncmp(inst->method, "Benders", 7) == 0){									// Benders' Loop method with time limit and patching heuristic
		double lb = 0;
		double ub = INFINITY;
		inst->succ = (int*)malloc(inst->nnodes * sizeof(int));
		inst->comp = (int*)malloc(inst->nnodes * sizeof(int));
		do{
			// Set time limit for cplex execution
			CPXsetdblparam(env, CPX_PARAM_TILIM, (inst->timelimit - (second()-inst->t_start)));
			
			// Set cutoff
			CPXsetdblparam(env, CPX_PARAM_CUTUP, ub);

			error = CPXmipopt(env,lp);
			if ( error ){
				printf("CPX error code %d\n", error);
				print_error("CPXmipopt() error"); 
			}
			double objval_p;
			CPXgetbestobjval(env, lp, &objval_p);
			lb = dmax(lb, objval_p);
			ncols = CPXgetnumcols(env, lp);
			xstar = (double *) calloc(ncols, sizeof(double));
				
			if ( CPXgetx(env, lp, xstar, 0, ncols-1) )
			{
				print_error("CPXgetx() error");
			}

			inst->ncols = ncols; // ncols corresponds to the number of variables in the model

			if (VERBOSE >= 50) {// print our solution				
				for (int i = 0; i < inst->nnodes; i++)
				{
					for (int j = i + 1; j < inst->nnodes; j++)
					{
						if (xstar[xpos(i, j, inst)] > 0.5) printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1); // x is in the solution if equal to 1

					}//for
				}//for
			}//if
			int* succ = (int*)malloc(inst->nnodes * sizeof(int)); // array with successors
			int* comp = (int*)malloc(inst->nnodes * sizeof(int)); // array with component numbers for each node
			int ncomp;											  //number of separate subtours
			build_sol(xstar, inst, succ, comp, &ncomp);

			if(ncomp > 1){
				updateModel(inst, env, lp, ncomp, comp);
				patchingHeuristicUpdate(env, lp, succ, comp, inst->nnodes, &ncomp, inst, NULL);
				double cost = calculate_succ_cost(succ, inst);
				if(cost < ub){
					ub = cost;
					inst->ncomp = ncomp;
					for (int k = 0; k < inst->nnodes; k++){
						inst->succ[k] = succ[k];
						inst->comp[k] = comp[k];
						}
				}
			}
			else{
				double cost = calculate_succ_cost(succ, inst);
				if(cost < ub){
					ub = cost;
					for (int k = 0; k < inst->nnodes; k++){
						inst->succ[k] = succ[k];
						inst->comp[k] = comp[k];
					}
					inst->ncomp = ncomp;
				}
				if(VERBOSE >= 30)
				{
					printf("Optimal\n");
					printf("LowerBound: %f; UpperBound: %f\n", lb, ub);
				}
				break;
			}
			if(VERBOSE >= 30)
			{
				printf("LowerBound: %f; UpperBound: %f\n", lb, ub);
			}
			free(succ);
			free(comp);
		}while((lb < (1-XSMALL)*ub) && ((second() - inst->t_start) < inst->timelimit));

	}else if(strncmp(inst->method, "Callback", 8) == 0){									// Callback
		int ncols = inst->ncols;

		//inizialization
		xstar = (double*)calloc(ncols, sizeof(double));	// here is my solution
		
		inst->succ = (int*)calloc(inst->nnodes, sizeof(int));
		inst->comp = (int*)calloc(inst->nnodes, sizeof(int));
		int ncomp;													// number of separate subtours

		if(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, my_callback, inst)) print_error("CPXcallbacksetfunc() error");
		
		extra_mileage(&inst->succ, 1, 1, 0, inst);
		
		//Converting a TSP solution succ[0..nnodes-1] to a CPLEX solution xheu[0..ncols-1]
		double* xheu = (double*)calloc(ncols, sizeof(double));
		for (int i = 0; i < ncols; i++) xheu[i] = 0;								//initialize xheu to all zeros
		for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, inst->succ[i], inst)] = 1;	//set solution edges to one
		int* ind = (int*)calloc(ncols, sizeof(int));
		for (int i = 0; i < ncols; i++) ind[i] = i; //indices
		int effortlevel = CPX_MIPSTART_NOCHECK;
		int beg = 0;

		//provide to cplex an incumbent solution
		if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");
		free(ind);
		free(xheu);
		error = CPXmipopt(env,lp);
		if ( error ) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}
		if ( CPXgetx(env, lp, xstar, 0, ncols-1) )
		{
			print_error("CPXgetx() error");
		}
		else {
			build_sol(xstar, inst, inst->succ, inst->comp, &ncomp);
			
			inst->ncomp = ncomp;
		}//if-else
	}else if(strncmp(inst->method, "Hard", 4) == 0){									// Math-Heuristic Hard-Fixing
		
		xstar = (double*)calloc(ncols, sizeof(double));
		extra_mileage(&inst->succ, 1, 1, 0, inst);
		double best_cost = INFINITY;
		ncols = CPXgetnumcols(env, lp);

		int* succ = (int*)malloc(inst->nnodes * sizeof(int));		// array with successors
		int* comp = (int*)malloc(inst->nnodes * sizeof(int));		// array with component numbers for each node
		inst->comp = (int*)malloc(inst->nnodes * sizeof(int));
		int ncomp;													// number of separate subtours
		double* xheu = (double*)calloc(ncols, sizeof(double));
		for (int i = 0; i < ncols; i++) xheu[i] = 0;								//initialize xheu to all zeros
		for (int i = 0; i < inst->nnodes; i++)
		{
			xheu[xpos(i, inst->succ[i], inst)] = 1;	//set solution edges to one
		}
		int* ind = (int*)calloc(ncols, sizeof(int));
		for (int i = 0; i < ncols; i++) ind[i] = i;
		int effortlevel = CPX_MIPSTART_NOCHECK;
		int beg = 0;
		double prob = inst->p;
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, my_callback, inst);
		
		if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");	//shifted from inside the loop
		do{
			freeAllLowerBound(env, lp, inst);
			for(int k=0; k<inst->nnodes; k++){
				if(random01() < prob && xheu[k] > 0.5){
					fixLowerBound(env, lp, k, inst->succ[k], inst);
				}
			}
			CPXsetdblparam(env, CPXPARAM_TimeLimit, (inst->timelimit - (second()-inst->t_start)));
			error = CPXmipopt(env,lp);
			if ( error ) 
			{
				printf("CPX error code %d\n", error);
				print_error("CPXmipopt() error"); 
			}
			if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("CPXgetx() error");	
			else{
				build_sol(xheu, inst, succ, comp, &ncomp);
				int current_cost = calculate_succ_cost(inst->succ, inst);
				if (current_cost < best_cost) {
					best_cost = current_cost;
					for (int k = 0; k < inst->nnodes; k++)
					{
						inst->succ[k] = succ[k];
						inst->comp[k] = comp[k];
					}
			}
			}
			calculateComponents(&inst->succ, &inst->comp, &inst->ncomp, xstar, inst);
		}while ((second()-inst->t_start) < inst->timelimit);
		inst->best_sol = &best_cost;
		free(xheu);
		free(ind);
		free(succ);
		free(comp);
	}else if(strncmp(inst->method, "Local", 5) == 0){									// Local Branching
		// Hyperparameter to set
		int k = 20;
		xstar = (double*)calloc(ncols, sizeof(double));

		// Find a starting solution with an heuristic:
		greedy(&inst->succ, 1, 1, 0, inst);
		two_opt(inst->succ, INFINITY, inst);
		double *xh = (double *) calloc(ncols, sizeof(double));
		for(int i=0; i<inst->nnodes; i++){
			if (i == inst->succ[i]) print_error(" error in succ list generated by the first heuristic algorithm");
			xh[xpos(i, inst->succ[i], inst)] = 1;
		}

		inst->ncols = CPXgetnumcols(env, lp);
		int ncols = inst->ncols;
		int n = inst->nnodes;
		int* succ = (int*)calloc(inst->nnodes, sizeof(int));		// array with successors
		int* comp = (int*)calloc(inst->nnodes, sizeof(int));		// array with component numbers for each node
		inst->comp = (int*)calloc(inst->nnodes, sizeof(int));
		int ncomp;		
		int* ind = (int*)calloc(ncols, sizeof(int));
		for (int i = 0; i < ncols; i++) ind[i] = i;
		int effortlevel = CPX_MIPSTART_NOCHECK;
		int beg = 0;

		if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xh, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");
		
		// Add this starting solution to the model
		localBranchingAddRow(env, lp, xh, k, inst);

		double step_time = inst->timelimit/4;

		// Set callback function
		if(CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, my_callback, inst)) print_error("CPXcallbacksetfunc() error");
		double best_cost = INFINITY;
		

		do{
			double remaining_time = step_time;
			if((inst->timelimit - (second()-inst->t_start))<step_time){
				remaining_time = (inst->timelimit - (second()-inst->t_start));
			}
			CPXsetdblparam(env, CPXPARAM_TimeLimit, remaining_time);
			if(VERBOSE >= 50){
				printf("cplex max time for this iteration: %f\n", remaining_time);
			}

			error = CPXmipopt(env,lp);
			if(error){
				printf("CPX error code %d\n", error);
				print_error("CPXmipopt() error");
			}
			if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("CPXgetx() error");
			calculateComponents(&inst->succ, &inst->comp, &inst->ncomp, xstar, inst);
			if(best_cost < calculate_succ_cost(inst->succ, inst)){
				print_error(" error in local branching, prev_cost can't be lower than the new one");
			}
			if(best_cost == calculate_succ_cost(inst->succ, inst)){
				k = k + 10;
				if(k > (inst->nnodes-1)/2) k = (inst->nnodes-1)/2;
			}
			best_cost = calculate_succ_cost(inst->succ, inst);
        	if(VERBOSE >= 50){
				printf("CPLEX Total cost: %lf\n", best_cost);
			}
			// Delete last row containing local branching constraint
			int nrows = CPXgetnumrows(env, lp);
			CPXdelrows(env, lp, nrows-1, nrows-1);

			// Create the new local branching constraint from the new solution found
			localBranchingAddRow(env, lp, xstar, k, inst);
			
		}while((second()-inst->t_start) < inst->timelimit);

		free(xh);
		free(ind);
		free(succ);
		free(comp);
	}

    if(VERBOSE >= 60){
        printf("conn comp algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->comp[i]);
        printf("\n");
    }
    if(VERBOSE >= 60){
        printf("Best solution found with CPLEX algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->succ[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("CPLEX Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
    }
	free(xstar);
	
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0; // or an appropriate nonzero error code

}

/***************************************************************************************************************************/
int xpos(int i, int j, instance *inst)                                           
/***************************************************************************************************************************/
{ 
	if ( i == j ){
		printf("warning: %d\n", i);
		print_error(" i == j in xpos" );
	}
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}


/***************************************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = dist(i,j,inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if ( VERBOSE >= 100 ) CPXwriteprob(env, lp, "model.lp", NULL);   

}


void localBranchingAddRow(CPXENVptr env, CPXLPptr lp, const double *xstar, int k, instance *inst){
	int *index = (int *) calloc(inst->nnodes * inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));
	char sense = 'G';
	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	int nnz = 0;
	double rhs = inst->nnodes - k;
	for(int i=0; i<inst->nnodes; i++){
		for(int j=i+1; j<inst->nnodes; j++){
			if(xstar[xpos(i,j,inst)] > 0.5){
				
				index[nnz] = xpos(i,j,inst);
				value[nnz] = 1.0;
				nnz ++;
			}
		}
	}
	int izero = 0;
	sprintf(cname[0], "localBranching");
	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");

}

// Slide 10; 18-May-2022 - MATH-HEURISTIC
void fixLowerBound(CPXENVptr env, CPXLPptr lp, int k, int succ_k, instance *inst){

	char lu = 'L';

	if(k == succ_k) print_error(" error in function fixLowerBound(...), k and successor of k(succ_k) can't be equal");
	int min, max;
	if(k>succ_k){
		min = succ_k;
		max = k;
	}else{
		min = k;
		max = succ_k;
	}

	int idx = xpos(min, max, inst);
	double lb = 1.0;
	if ( CPXchgbds(env, lp, 1, &idx, &lu, &lb) ) print_error(" wrong CPXnewcols on x var.s");

}

// Slide 10; 18-May-2022 - MATH-HEURISTIC
void freeLowerBound(CPXENVptr env, CPXLPptr lp, int min_k, int max_k, instance *inst){

	char lu = 'L';

	if(min_k >= max_k) print_error(" error in function fixLowerBound(...), min_k and can't be greater or equal than max_k");
	int idx = xpos(min_k, max_k, inst);
	double lb = 0.0;
	if ( CPXchgbds(env, lp, 1, &idx, &lu, &lb) ) print_error(" wrong CPXnewcols on x var.s");
}

// Slide 10; 18-May-2022 - MATH-HEURISTIC
void freeAllLowerBound(CPXENVptr env, CPXLPptr lp, instance *inst){
	for(int i=0; i<inst->nnodes; i++){
		for(int j=i+1; j<inst->nnodes; j++){
			freeLowerBound(env, lp, i, j, inst);
		}
	}
}

// Slide --; 4-May-2022 - CALLBACK
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) 
{ 
	
	instance* inst = (instance*) userhandle;  
	double* xstar = (double*) malloc(inst->ncols * sizeof(double));  
	double objval = CPX_INFBOUND; 
	if ( CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols-1, &objval) ) print_error("CPXcallbackgetcandidatepoint error");
	
	// get some random information at the node (as an example for the students)
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread); 
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode); 
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent); 

	int nnz = 0; 
	int ncomp = 10;
	int* comp; // array with component numbers for each node
	int *a;
	
	calculateComponents(&a, &comp, &ncomp, xstar, inst);

	if(ncomp>1){
		// Add SEC's
		for(int k=0; k<ncomp; k++){
			addOneCompSec(NULL, NULL, comp, k, context, inst);
		}
		// Patching heuristic
		if(inst->patch){
			free(comp);
			calculateComponents(&a, &comp, &ncomp, xstar, inst);
			patchingHeuristicUpdate(NULL, NULL, a, comp, inst->nnodes, &ncomp, inst, context);
		}
	}

	double objheu = calculate_succ_cost(a, inst);
	incumbent = CPX_INFBOUND;
	CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);

	//Post the solution to CPLEX if with no subtours
	if(inst->post && ncomp==1 && objheu<(1-XSMALL)*incumbent){
		double* xheu = (double*)malloc(inst->ncols * sizeof(double));

		//Converting a "successor" TSP solution succ[0..nnodes-1] to a CPLEX solution xheu[0..ncols-1]
		for (int i = 0; i < inst->ncols; i++) xheu[i] = 0.0;							//initialize xheu to all zeros
		for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, a[i], inst)] = 1.0;	//set solution edges to one
		int* ind = (int*)malloc(inst->ncols * sizeof(int));
		for (int j = 0; j < inst->ncols; j++) ind[j] = j;

		two_opt(a, inst->ref_timelimit, inst);
		if (CPXcallbackpostheursoln(context, inst->ncols, ind, xheu, calculate_succ_cost(a, inst), CPXCALLBACKSOLUTION_NOCHECK)) print_error("CPXcallbackpostheursoln() error");

		//free memory
		free(ind);
		free(xheu);
	}
	free(a);
	free(comp);
	free(xstar); 
	return 0; 
}


// Slide 06; 13-Apr-2022
void updateModel(instance *inst, CPXENVptr env, CPXLPptr lp, int ncomp, int* comp){

	// Adding SEC's to CPLEX's model
	for(int k = 0; k < ncomp; k++){
		addOneCompSec(env, lp, comp, k, NULL, inst);
	}
}

// Slide 06; 13-Apr-2022
void addOneCompSec(CPXENVptr env, CPXLPptr lp, int *comp, int comp_number, CPXCALLBACKCONTEXTptr context, instance *inst){

	int *index = (int *) calloc(inst->nnodes * inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));
	char sense = 'L';
	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	int nnz = 0;
	double rhs = -1.0;
	// For each vertex with his component equals to k
	for(int i=0; i<inst->nnodes; i++){
		if(comp[i] != comp_number) continue;
		rhs ++;
		for(int j=i+1; j<inst->nnodes; j++){
			if(comp[j] != comp_number) continue;

			index[nnz] = xpos(i,j,inst);
			value[nnz] = 1.0;
			nnz ++;
		}
	}
	int izero = 0;
	sprintf(cname[0], "sec");
	if(context == NULL){
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	}else{
		if ( CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value) ) 
			print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
	}
	free(index);
	free(value);
	
}

// Slide 06; 13-Apr-2022
void calculateComponents(int **succ, int **comp, int *ncomp, const double *xstar, instance *inst){

	(*succ) = (int *) calloc(inst->nnodes, sizeof(int));
	(*comp) = (int *) calloc(inst->nnodes, sizeof(int));
	for(int i=0; i<inst->nnodes; i++) (*succ)[i] = -1;
	for(int i=0; i<inst->nnodes; i++) (*comp)[i] = -1;

	int num_comp = -1;

    for (int i = 0; i < inst->nnodes; i++){
        if ((*comp)[i] >= 0)
            continue;
        num_comp++;
        int next = i;
        int done = 0;
        while (!done){
            (*comp)[next] = num_comp;
            done = 1;
            for (int j = i+1; j < inst->nnodes; j++){
                if (next == j || (*comp)[j] >= 0)
                    continue;
                if (xstar[xpos(next, j, inst)] >= 0.5){
                    (*succ)[next] = j;
                    next = j;
                    done = 0;
                    break;
                }
            }
        }
        (*succ)[next] = i;
    }
	*ncomp = num_comp+1;
}

// Slide 07; 22-Apr-2022
void patchingHeuristicUpdate(CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int nnodes, int *ncomp, instance *inst, CPXCALLBACKCONTEXTptr context){
	
	int oldcomp = *ncomp;
	for(int k=0; k<oldcomp-1; k++){
		int a, next_a;
		int b, next_b;
		double delta = INFINITY;
		// For each node i check for another node j with different comp s.t. delta is the minimum possible
		for(int i=0; i<nnodes; i++){
			for(int j=0; j<nnodes; j++){
				if(comp[i] >= comp[j]) continue; // if i == j they have the same comp[].
				if(succ[i] == j) continue;
				if(succ[j] == i) continue;
				double temp_delta = (cost(i, succ[j], inst) + cost(j, succ[i], inst))-(cost(i, succ[i], inst) + cost(j, succ[j], inst));//new_solution - old_solution
				if(temp_delta<delta){
					a = i;
					next_a = succ[i];
					b = j;
					next_b = succ[j];
					delta = temp_delta;
				}
			}
		}
		// Merge 2 components and keep the minimum which is on 'a' since "if(comp[i] >= comp[j]) continue;" -> comp[a] < comp[b]
		succ[a] = next_b;
		succ[b] = next_a;
		// From next_b to b change comp[]
		int prev = next_b;
		int next = succ[prev];
		comp[prev] = comp[a];
		int counter = 0;
		while(succ[prev] != b){
			comp[next] = comp[a];
			prev = next;
			next = succ[next];
			if(counter >= inst->nnodes){
				print_error(" loop on patchingHeuristic not working");
			}
			counter ++;
		}
		comp[next] = comp[a];
		*ncomp = *ncomp-1;
		if(k==(oldcomp-2)) break;
		
		// Add SEC's
		addOneCompSec(env, lp, comp, comp[a], context, inst);
	}
}



#define DEBUG    // comment out to avoid debugging 
#define EPS 1e-5

/*********************************************************************************************************************************/
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{   
#ifdef DEBUG
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ((fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0) > EPS))
			{
				print_error(" wrong xstar in build_sol()");
			}
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);
#endif

	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp-1;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}


	
	
	





 






