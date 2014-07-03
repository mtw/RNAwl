/*  Last changed Time-stamp: <2014-07-03 09:43:27 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

move_str get_random_move_pt(const char *,const short int*);
void apply_move_pt(short int *,const move_str);

#endif
