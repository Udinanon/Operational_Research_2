#include "utility.h"

static enum LOG_LEVEL log_level = DEBUG;
static FILE* log_file = NULL;

// https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html
void logger(enum LOG_LEVEL level, const char *fmt, ...) {
    if (level > log_level){
        return;
    }
    va_list ap;
    time_t t;
    char datestr[51];

    /* determine if we just go to std error */
    FILE* file = (log_file == NULL) ? stderr : log_file;

    /* datetime & pid formatting */
    t = time(NULL);
    tzset();
    strftime(datestr, sizeof(datestr) - 1, "%a %b %d %T %Z %Y", localtime(&t));
    fprintf(file, "%s [%d]: ", datestr, getpid());

    /* draw out the vararg format */
    va_start(ap, fmt);
    vfprintf(file, fmt, ap);
    va_end(ap);

    /* bump to the next line */
    fprintf(file, "\n");
}

void set_logger_level(enum LOG_LEVEL level){
    log_level = level;
}

int set_logger_file(char* filename){
    log_file = fopen(filename, "r+");
    if (log_file == NULL) {  // check for file access,
        logger(ERROR, "Failed opening log file %s", filename);
        perror("ERROR: ");
        return -1;
    }
    return 0;
}

double get_random(){
    return (double)((double)rand() / RAND_MAX);
}