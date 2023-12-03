#include "tsp.h"


int two_opt(int *succ, double t_limit, instance *inst){
    double t_start = second();
    double delta;
    int temp_i, temp_j;
    
    do{
        delta = INFINITY;
        for(int i=0; i<inst->nnodes-1; i++){
            for(int j=i+1; j<inst->nnodes; j++){
                double temp_delta = (cost(i,j,inst) + cost(inst->succ[i], inst->succ[j], inst)) - (cost(i, inst->succ[i], inst) + cost(j, inst->succ[j], inst));
                if(temp_delta<delta){
                    delta = temp_delta;
                    temp_i = i;
                    temp_j = j;
                }
            }
        }

        if(VERBOSE >= 80){
            printf("temp_i: %d, temp_j: %d\n",temp_i, temp_j);
            printf("delta: %f\n",delta);
        }
        if(delta<0){
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
        }

    }while((delta<0) && ((second()-t_start)<t_limit));

    if(VERBOSE >= 40){
        printf("Best solution found with 2-OPT algorithm:\n");
        for(int i=0; i<inst->nnodes; i++) printf("|%d=%d| ",i, inst->succ[i]);
        printf("\n");
    }
    if(VERBOSE >= 30){
        printf("2-OPT Total cost: %lf\n", calculate_succ_cost(inst->succ, inst));
    }
    return 0;
}