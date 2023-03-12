#include "tsp_parser.h"

static int compare_string(const char* s1, const char* s2);

void parse_file(const char* filename, TSP_data* data) {
    FILE* fp;  // open file
    logger(INFO, "Opening file %s", filename);
    if ((fp = fopen(filename, "r+")) == NULL) {  // check for file access,
        logger(ERROR, "ERROR READING FILE %s", filename);
        perror("ERROR: ");
        exit(EXIT_FAILURE);
    }
    logger(INFO, "Opened file");
    // allocate variables to read filelines
    char* line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = -1;
    int line_size = 0;

    int data_section = 0;
    char* label = NULL;
    int data_counter = 0;

    // read lines in order
    while ((line_size >= 0)) {
        line_count++;
        line_size = getline(&line_buf, &line_buf_size, fp);
        if (line_size < 0 ) break;
        logger(ALL, "Reached line #%d: %s", line_count, line_buf);
        if (data_section == 0){       
            label = strtok(line_buf, ":");

            if (compare_string(label, "NAME")){
                logger(INFO, "PROBLEM NAME: %s", strtok(NULL, ":"));
            }
            if (compare_string(label, "COMMENT")) {
                logger(INFO, "PROBLEM COMMENT: %s", strtok(NULL, ":"));
            }
            if (compare_string(label, "TYPE")) {
                logger(INFO, "PROBLEM TYPE: %s", strtok(NULL, ":"));
            }
            if (compare_string(label, "DIMENSION")) {
                int dimensions = strtol(strtok(NULL, ":"), NULL, 10);
                logger(INFO, "PROBLEM DIMENSION: %d\n", dimensions);
                data->n_dimensions = dimensions;
                allocate_points(data);
            }
            if (compare_string(label, "EDGE_WEIGHT_TYPE")) {
                logger(DEBUG, "PROBLEM EDGE_WEIGHT_TYPE: %s", strtok(NULL, ":"));
            }
            if (compare_string(label, "NODE_COORD_SECTION")) {
                logger(INFO, "REACHED NODE_COORD_SECTION\n");
                if (data->n_dimensions == -1){
                    logger(FATAL, "Reached Coords section, never read dimensions. File corrupt");
                    exit(EXIT_FAILURE);
                }
                data_section = 1;
                continue;
            }
        }
        if (data_section == 1){
            char* section = strtok(line_buf, " ");
            if (compare_string(section, "EOF")){
                logger(INFO, "Reached EOF mark. Counter %d dimension %d", data_counter, data->n_dimensions);
                break;
            }
            long index = strtol(section, NULL, 10);
            double x = strtod(strtok(NULL, " "), NULL);
            double y = strtod(strtok(NULL, " "), NULL);
            if (index - data_counter != 1){
                logger(WARN, "MISALIGNMENT IN DATA: INDEX %d COUNTER %d", index, data_counter);
            }
            //data->points[data_counter] = malloc(sizeof(Point));
            Point* point = &data->points[data_counter];
            point->index = index;
            point->x = x;
            point->y = y;
            data_counter++;
            logger(ALL, "STORED DATA: %d %f %f", point->index, point->x, point->y);
        }
    }
    fclose(fp);
    free(line_buf);
}

static int compare_string(const char* s1,const  char* s2){
    return (strncmp(s1, s2, strlen(s2)) == 0);
}