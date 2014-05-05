#ifndef GLOBALS_H
#define GLOBALS_H

//#include "config.h"
#include <gsl/gsl_histogram.h>

/* variables */
float mfe;

/* function pointers */
void  (*pre_process_model)(void);
void  (*post_process_model)(void);
//float (*energie)(const char *, const char *);
void  (*initialize_model)(char *);
//void  (*move_it)(char *);

/* functions */
void process_commandline (int argc, char *argv[]);
void wanglandau(void);
void wanglandau_free_memory(void);
void sighandler (int);


/* RNA-related */
void initialize_RNA(char *);
void pre_process_RNA(void);
void post_process_RNA(void);
//void subopt_first_bin_RNA(double);



#endif
