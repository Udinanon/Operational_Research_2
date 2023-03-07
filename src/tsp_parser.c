#include "tsp_parser.h"

static int compare_string(const char* s1, const char* s2);

TSP_data parse_file(char* filename) {

    FILE* fp;  // open file
    if ((fp = fopen(filename, "r+")) == NULL) {  // check for file access,
        printf("ERROR READING FILE %s", filename);
        perror("ERROR: ");
    }
    printf("Acessed file\n");
    // allocate variables to read filelines
    char* line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = 0;
    int line_size;

    TSP_data tsp_data;
    int data_section = 0;
    char* label = NULL;
    int data_counter = 0;

    // read lines in order
    line_size = getline(&line_buf, &line_buf_size, fp);
    while ((line_size >= 0)) {
        if (data_section == 0){       
            label = strtok(line_buf, ":");

            if (compare_string(label, "NAME")){
                printf("PROBLEM NAME: %s\n", strtok(NULL, ":"));
            }
            if (compare_string(label, "COMMENT")) {
                printf("PROBLEM COMMENT: %s\n", strtok(NULL, ":"));
            }
            if (compare_string(label, "TYPE")) {
                printf("PROBLEM TYPE: %s\n", strtok(NULL, ":"));
            }
            if (compare_string(label, "DIMENSION")) {
                int dimensions = strtol(strtok(NULL, ":"), NULL, 10);
                printf("PROBLEM DIMENSION: %d\n", dimensions);
                tsp_data.n_dimensions = dimensions;
                tsp_data.points = malloc(sizeof(Point) * dimensions);
            }
            if (compare_string(label, "EDGE_WEIGHT_TYPE")) {
                printf("PROBLEM EDGE_WEIGHT_TYPE: %s\n", strtok(NULL, ":"));
            }
            if (compare_string(label, "NODE_COORD_SECTION")) {
                printf("REACHED NODE_COORD_SECTION\n");
                data_section = 1;
            }
        }
        if (data_section == 1){
            long index = strtol(strtok(&line_buf, " "), NULL, 10);
            double x = strtod(strtok(NULL, " "), NULL);
            double y = strtod(strtok(NULL, " "), NULL);
            if (index - data_counter != 1){
                printf("MISALIGNMENT IN DATA: INDEX %d COUNTER %d\n", index, data_counter);
            }
            Point* point = tsp_data.points[data_counter];
            point->index = index;
            point->x = x;
            point->y = y;
            data_counter++;
        }



        // prepare next line
        line_count++;
        line_size = getline(&line_buf, &line_buf_size, fp);
    }
    free(fp);
    free(line_buf);
    return tsp_data;
}

static int compare_string(const char* s1,const  char* s2){
    return (strncmp(s1, s2, strlen(s2)) == 0);
}