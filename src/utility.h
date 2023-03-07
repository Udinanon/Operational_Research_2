#ifndef UTILITY_H

#define UTILITY_H

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

enum LOG_LEVEL{FATAL, ERROR, WARN, INFO, DEBUG, ALL};

/** Logs a line of text */
void logger(enum LOG_LEVEL level, const char *fmt, ...);
void set_logger_level(enum LOG_LEVEL level);
int set_logger_file(char* filename);

#endif