/*
  globals.h : global definitions for Wang-Landau sampling
  Last changed Time-stamp: <2014-07-02 16:54:48 mtw>
*/

#ifndef GLOBALS_H
#define GLOBALS_H

#include "config.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

/* variables/arrays */
float mfe;
gsl_histogram *h;    /* histogram of energies seen in current iteration */
gsl_histogram *s;    /* true DOS of the lowest energy range (if
			available); required for normalization, which
			is computed based on the lowest-energy bins */
double *dos;         /* scaled and normalized DOS */

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
void dealloc_gengetopt (void);
#endif
