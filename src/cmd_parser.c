#include "cmd_parser.h"

static int parse_opt(int key, char* arg, struct argp_state* state);
static void prepare_params(Parameters*);

Parameters parse(int counter, char** args){
    Parameters params;
    prepare_params(&params);

    const char* argp_program_version = "TSPSolver 0.0.0";
    const char* argp_program_bug_address = "<martino.trapanotto@studenti.unipd.it>";

    struct argp_option options[] = {
        {"tsp-file", 'f', "TSP_FILENAME", 0, "TSPlib file to load"},
        {"verbosity", 'v', "LEVEL", OPTION_ARG_OPTIONAL, "Level of verbosity in logging. Default to ALL"},
        {"seed", 's', "SEED", OPTION_ARG_OPTIONAL, "Randomness seed. Default is random"},
        {"n-thread", 'n', "THREADS", OPTION_ARG_OPTIONAL, "Number of threads. Default is 1"},
        {0}
    };
    struct argp argp = {options, parse_opt};
    argp_parse(&argp, counter, args, 0, 0, &params);
    return params;
}

static int parse_opt(int key, char* arg, struct argp_state* state) {
    struct Parameters* params = state->input;
    switch (key) {
        case 'f':
            params->filename = arg;
            break;
        case 'v':
            params->verbosity = (enum LOG_LEVEL)strtol(arg, NULL, 10);
            break;
        case 's':
            params->seed = strtol(arg, NULL, 10);
            break;
        case 'n':
            params->n_threads = strtol(arg, NULL, 10);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static void prepare_params(Parameters* params){
    params->n_threads = 1;
    params->verbosity = ALL;
}

static int check_parms(Parameters* params){
    if (params->filename == NULL){
        return -1;
    }
    return 0;
}