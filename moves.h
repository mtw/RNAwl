/*  Last changed Time-stamp: <2014-05-09 23:43:06 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

move_str get_random_move_pt(const char *,const short int*,int);
void apply_move_pt(short int *,const move_str);

#endif
