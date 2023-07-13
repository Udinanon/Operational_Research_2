#include "tsp.h"

void parse_command_line(int argc, char** argv, instance *inst) { 
	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	inst->model_type = 0;
	inst->method = 0;
	inst->p = 1;
	inst->len_rcl = 1;
	inst->refinement = -1;
	inst->meta = -1;
	inst->kick = -1;
	inst->population = -1;
	strcpy(inst->input_file, "NULL");
	// strcpy(inst->output_file, "NULL");
	inst->randomseed = 0; 
	inst->num_threads = 0;
	inst->timelimit = INFINITY; 
	inst->nrand = 0;
	inst->integer_costs = 0;
	inst->plot = 1;
	inst->post = 0;

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		// if ( strcmp(argv[i],"-out_file") == 0 ) { strcpy(inst->output_file,argv[++i]); continue; } 		// output file
		// if ( strcmp(argv[i],"-of") == 0 ) { strcpy(inst->output_file,argv[++i]); continue; } 			// output file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }				// total time limit
		if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-method") == 0 ) { inst->method = atoi(argv[++i]); continue; } 			// method for cplex
		if ( strcmp(argv[i],"-p") == 0 ) { inst->p = atof(argv[++i]); continue; } 						// probability param
		if ( strcmp(argv[i],"-len_rcl") == 0 ) { inst->len_rcl = atoi(argv[++i]); continue; } 			// length of RCL param
		if ( strcmp(argv[i],"-refinement") == 0 ) { inst->refinement = atoi(argv[++i]); continue; } 	// refinement type
		if ( strcmp(argv[i],"-meta") == 0 ) { inst->meta = atoi(argv[++i]); continue; } 				// meta type
		if ( strcmp(argv[i],"-kick") == 0 ) { inst->kick = atoi(argv[++i]); continue; } 				// kick param
		if ( strcmp(argv[i],"-population") == 0 ) { inst->population = atoi(argv[++i]); continue; } 				// kick param
		if ( strcmp(argv[i],"-patch") == 0 ) { inst->patch = atoi(argv[++i]); continue; }				// patchheuristic
		if ( strcmp(argv[i],"-post") == 0 ) { inst->post = atoi(argv[++i]); continue; }					// postheuristic
		if ( strcmp(argv[i],"-n") == 0 ) { inst->nrand = abs(atoi(argv[++i])); continue; } 				// number of random points
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-plot") == 0 ) { inst->plot = atoi(argv[++i]); continue; } 		// plot or not
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters (vers. 16-may-2015) --------------------------------------------------\n");
		printf("-file || -f %s\n", inst->input_file); 
		printf("-time_limit || -tl %lf\n", inst->timelimit); 
		printf("-model_type || -model %d\n", inst->model_type); 
		printf("-method %d\n", inst->method); 
		printf("-p %f\n", inst->p); 
		printf("-len_rcl %d\n", inst->len_rcl); 
		printf("-refinement %d\n", inst->refinement); 
		printf("-meta %d\n", inst->meta); 
		printf("-kick %d\n", inst->kick); 
		printf("-population %d\n", inst->population); 
		printf("-n %d\n", inst->nrand); 
		printf("-seed %d\n", inst->randomseed); 
		printf("-threads %d\n", inst->num_threads); 
		printf("-int %d\n", inst->integer_costs); 
		printf("-patch %d\n", inst->patch);
		printf("-post %d\n", inst->post);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if ( help ) exit(1);

}


void check_feasibility(instance *inst){	//checker

	int prev = 0;
	int succ;
	int *unv_nodes = (int *) calloc(inst->nnodes, sizeof(int));
	int len_unv_nodes = inst->nnodes;
	succ = inst->succ[prev];
	while(len_unv_nodes>0){

		if(unv_nodes[prev] == 1){
			if(len_unv_nodes==0 && prev == 0)
			{
				printf("Printing last node");
			}
			else{
				print_error("Found inner cycle, not feasible solution!");
			}
		}
		unv_nodes[prev] = 1;
		prev = inst->succ[prev];
		succ = inst->succ[prev];
		len_unv_nodes --;
	}
}



void plot_best_sol(instance *inst){
	FILE *fout = fopen("../plot/bs2.dat", "w");
	if(fout == NULL) print_error(" can't write to output file");

	int prev = 0;
	int succ;
	int *unv_nodes = (int *) calloc(inst->nnodes, sizeof(int));
	int len_unv_nodes = inst->nnodes;

	fprintf(fout, "%lf %lf 0\r\n", inst->xcoord[prev], inst->ycoord[prev]);
	succ = inst->succ[prev];
	fprintf(fout, "%lf %lf 0\r\n", inst->xcoord[succ], inst->ycoord[succ]);
	fprintf(fout, "\r\n");
	while(len_unv_nodes>0){

		if(unv_nodes[prev] == 1){
			if(len_unv_nodes==0 && prev == 0)
			{
				printf("Printing last node");
			}
			else{
				print_error("Found inner cycle, not feasible solution!");
				if(len_unv_nodes==0 && prev != 0)
				{
					prev = 0;
				}
				else{
					prev = -1;
					for(int i=0; (i<inst->nnodes) && (prev == -1); i++){
						if(unv_nodes[i] == 0) prev = i;
					}
				}
			}
			fprintf(fout, "%lf %lf 1\r\n", inst->xcoord[prev], inst->ycoord[prev]);
			succ = inst->succ[prev];
			fprintf(fout, "%lf %lf 1\r\n", inst->xcoord[succ], inst->ycoord[succ]);
			fprintf(fout, "\r\n");

		}
		unv_nodes[prev] = 1;
		prev = inst->succ[prev];
		fprintf(fout, "%lf %lf 1\r\n", inst->xcoord[prev], inst->ycoord[prev]);
		succ = inst->succ[prev];
		fprintf(fout, "%lf %lf 1\r\n", inst->xcoord[succ], inst->ycoord[succ]);
		fprintf(fout, "\r\n");
		len_unv_nodes --;
	}
	fclose(fout);

	system("gnuplot ../plot/settings_bs2.txt");
}


void plot_best_sol_old(instance *inst){
	FILE *fout = fopen("../plot/bs.dat", "w");
	if(fout == NULL) print_error(" can't write to output file");

	int prev = 0;
	int end = 0;
	fprintf(fout, "%lf %lf\r\n", inst->xcoord[prev], inst->ycoord[prev]);
	while(!end){
		prev = inst->succ[prev];
		//printf("Writing node: %d\n", prev);
		fprintf(fout, "%lf %lf\r\n", inst->xcoord[prev], inst->ycoord[prev]);
		if(prev == 0) end = 1;
	}
	fclose(fout);

	system("gnuplot ../plot/settings_bs.txt");
}


void plot_points(instance *inst){
	FILE *fout = fopen("../plot/data.dat", "w");
	if(fout == NULL) print_error(" can't write the output file");

	for(int i=0; i<inst->nnodes; i++){
		fprintf(fout, "%lf %lf\r\n", inst->xcoord[i], inst->ycoord[i]);
	}
	fclose(fout);

	system("gnuplot ../plot/settings.txt");
}


void read_input(instance *inst) // simplified TSP parser, not all SECTIONs detected  
{
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");
	
	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		else if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		else if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}
		

		else if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));    
			active_section = 0;  
			continue;
		}


		else if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{
			token1 = strtok(NULL, " :");
            if (strncmp(token1, "EUC_2D", 6) == 0)
                inst->weight_type = EUC_2D;
            if (strncmp(token1, "MAX_2D", 6) == 0)
                inst->weight_type = MAX_2D;
            if (strncmp(token1, "MAN_2D", 6) == 0)
                inst->weight_type = MAN_2D;
            if (strncmp(token1, "CEIL_2D", 7) == 0)
                inst->weight_type = CEIL_2D;
            if (strncmp(token1, "ATT", 3) == 0)
                inst->weight_type = ATT;
            if (strncmp(token1, "EXPLICIT", 8) == 0)
                print_error("#PARSE_INST Wrong edge weight type, this program resolve only 2D TSP case with coordinate type.");
			active_section = 0;
			continue;
		}            
		
		else if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}

		else if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}

		else{
		if ( strncmp(par_name, "DEMAND_SECTION", 14) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before DEMAND_SECTION section");
			active_section = 2;
			continue;
		}  
		
		if (strncmp(par_name, "EDGE_WEIGHT_SECTION", 19) == 0)
        {
			active_section = 2;
            continue;
        }
		}
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}

		if( active_section == 2){
			
			continue;
		}


		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");     
		    
	}                

	fclose(fin);    
	
}


void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }  

void delete_instance(instance* inst)
{
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->costs);
	free(inst->succ);
	free(inst->path);	
	free(inst->comp);
	free(inst->best_sol);
	free(inst);
	return;
}
