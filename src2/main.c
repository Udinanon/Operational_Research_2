#include "tsp.h"


void compute_costs(instance *inst);
int TSPopt(instance *inst);

int main(int argc, char **argv) 
{
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }


	double t1 = second(); 
	instance inst;

	parse_command_line(argc,argv, &inst);     
	
	// //printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);
	double *xc;
	double *yc;
	if(strcmp(inst.input_file, "NULL") != 0){
		read_input(&inst);
	}else if(inst.nrand != 0){
		generate_random_points(&inst.xcoord, &inst.ycoord, inst.nrand, inst.randomseed, 1000, 1000);
		inst.nnodes = inst.nrand;
	}

	inst.t_start = t1;
	compute_costs(&inst);

	if ( TSPopt(&inst) ) print_error(" error within TSPopt()");
	double t2 = second(); 
	if(inst.plot == 1){
    	plot_best_sol(&inst);
	}
	if ( VERBOSE >= 1 )   
	{
		printf("... TSP problem solved in %lf sec.s\n", t2-t1);  
	}
	
	free(xc);
	free(yc);
	delete_instance(&inst);
	return 0; 
}




