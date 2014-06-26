/* Last changed Time-stamp: <2014-06-26 17:10:11 mtw> */

#ifndef WL_OPTIONS_H
#define WL_OPTIONS_H

#include <stdio.h>

typedef struct _options {
  FILE *INFILE;        /* input file */
  char *basename;       /* base name of processed file */
  char *sequence;      /* sequence */
  char *structure;     /* start structure */
  int len;             /* sequence length */
  int bins;            /* # of equidistant bins in histogram */
  double ffinal;       /* modification parameter f */
  float flat;          /* flatness criterion */
  int steps;           /* wl steps before histogram is checked for flatness */
  float T;             /* fold temperature */
  float erange;        /* energy range for subopt */
  int norm;            /* # of normalization-bins */
  int verbose;         /* be verbose */
} options;

options wanglandau_opt;

#endif
