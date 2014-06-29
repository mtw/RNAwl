/*
  globals.h : global definitions for Wang-Landau sampling
  Last changed Time-stamp: <2014-06-29 00:49:30 mtw>
*/

#ifndef GLOBALS_H
#define GLOBALS_H

#include "config.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

/* variables/arrays */
float mfe;
gsl_histogram *h;     /* histogram */

/* function pointers */
void  (*pre_process_model)(void);
void  (*post_process_model)(void);
//float (*energie)(const char *, const char *);
void  (*initialize_model)(const char *);

/* functions */
void process_commandline (int argc, char *argv[]);
void wanglandau(void);
void wanglandau_free_memory(void);
void sighandler (int);

#endif
