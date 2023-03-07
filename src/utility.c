#include "utility.h"

static enum LOG_LEVEL log_level = DEBUG;

// https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html
void logger(FILE* file, enum LOG_LEVEL level, const char *fmt, ...) {
    if (level > log_level){
        return;
    }
    va_list ap;
    time_t t;
    char datestr[51];

    /* determine if we just go to std error */
    file = (file == NULL) ? stderr : file;

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

void setup_logger(enum LOG_LEVEL level){
    log_level = level;
}