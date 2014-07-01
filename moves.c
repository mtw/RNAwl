#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "globals.h"
#include "moves.h"
#include <ViennaRNA/move_set.h>

#define MINGAP 3
  
//int get_list(struct_en*, struct_en*);
static int construct_moves_new(const char*, const short*, int , move_str **);
inline int try_insert_seq2(const char*, int, int);
inline int compat(const char, const char);
void mtw_dump_pt(const short*);

/*
  compute a random move on a pair table
  returns move operations to be applied to pt in order to perform the move
 */
move_str
get_random_move_pt(const char *seq, const short int *pt,int verbose)
{
  move_str r,*mvs=NULL;
  int i,count;
  
  count = construct_moves_new((const char *)seq,pt,1,&mvs);
  
  if (verbose ==1){
    /*
    for (i = 0; i<count; i++) {  
      printf("%d %d\n", mvs[i].left, mvs[i].right);
    }
    */
    printf ("++ applying move: left %i right %i\n",mvs[0].left,mvs[0].right);
  }
  
  r.left  = mvs[0].left;
  r.right = mvs[0].right;
  free(mvs);
  return r;
}

/*
  apply move operation on a pair table
*/
void
apply_move_pt(short int *pt,
	      move_str m)
{
  if(m.left < 0){
    pt[(int)(fabs(m.left))] = 0;
    pt[(int)(fabs(m.right))] = 0;
  }
  else {
    pt[m.left] = m.right;
    pt[m.right] = m.left;
  }
  //print_str(stdout,pt);printf("\n");
}

static int
construct_moves_new(const char *seq,
		    const short *structure,
		    int permute,
		    move_str **array)
{
  /* generate all possible moves (less than n^2)*/
  int i;
  int size = 4;
  int count = 0;
  move_str *res = (move_str*) malloc(sizeof(move_str)*(size+1));
  
  for (i=1; i<=structure[0]; i++) {
    if (structure[i]!=0) {
      if (structure[i]<i) continue;
      count ++;
      // need to reallocate the array?
      if (count>size) {
      	size *= 2;
      	res = realloc(res, sizeof(move_str)*(size));
      }
      res[count-1].left = -i;
      res[count-1].right = -structure[i];
      //fprintf(stderr, "add  d(%d, %d)\n", i, structure[i]);
    } else {
      int j;
      for (j=i+1; j<=structure[0]; j++) {
        /* fprintf(stderr, "check (%d, %d)\n", i, j); */
        if (structure[j]==0) {
          if (try_insert_seq2(seq,i,j)) {
            count ++;
	    // need to reallocate the array?
	    if (count>size) {
	      size *= 2;
	      res = realloc(res, sizeof(move_str)*(size));
	    }
	    res[count-1].left = i;
	    res[count-1].right = j;
            //fprintf(stderr, "add  i(%d, %d)\n", i, j);
            continue;
          }
        } else if (structure[j]>j) { /*  '(' */
          j = structure[j];
        } else break;
      }
    }
  }
  
  res = realloc(res, sizeof(move_str)*(count));
  
  /* permute them */
  if (permute) {
    for (i=0; i<count; i++) {
      int rnd = rand(); 
      rnd = rnd % (count-i) + i;
      move_str mv;
      mv = res[i];
      res[i] = res[rnd];
      res[rnd] = mv;
    }
  }
  *array = res;
  return count;
}

/*  try insert base pair (i,j) */
inline int
try_insert_seq2(const char *seq,
	       int i,
	       int j)
{
  if (i<=0 || j<=0) return 0;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]));
}

/* compatible base pair?*/
inline int
compat(const char a,
       const char b)
{
  if (a=='A' && b=='U') return 1;
  if (a=='C' && b=='G') return 1;
  if (a=='G' && b=='U') return 1;
  if (a=='U' && b=='A') return 1;
  if (a=='G' && b=='C') return 1;
  if (a=='U' && b=='G') return 1;
  /* and with T's*/
  if (a=='A' && b=='T') return 1;
  if (a=='T' && b=='A') return 1;
  if (a=='G' && b=='T') return 1;
  if (a=='T' && b=='G') return 1;
  return 0;
}


void
mtw_dump_pt(const short *pairtable)
{
  int i;
  printf("> ");
  for (i=0;i<*pairtable;i++){
    printf("%i ",*(pairtable+i));
  }
  printf("\n");
}

