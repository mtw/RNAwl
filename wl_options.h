#ifndef WL_OPTIONS_H
#define WL_OPTIONS_H

#include <stdio.h>

typedef struct _options {
  FILE *INFILE;
  char *sequence;
  char *structure;
  int len;
  int bins;
  float erange;
} options;

options wanglandau_opt;

#endif
