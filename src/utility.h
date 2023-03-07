#ifndef UTILITY_H

#define UTILITY_H

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

enum LOG_LEVEL{FATAL, ERROR, WARN, INFO, DEBUG, ALL};

/** Logs a line of text */
void logger(FILE *file, enum LOG_LEVEL level, const char *fmt, ...);


#endif