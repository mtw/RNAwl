#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "moves.h"
#include <gsl/gsl_histogram.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/pair_mat.h>

#define MINGAP 3

typedef struct move_str {
  int left;
  int right;
} move_str;

typedef struct _options {
  FILE *INFILE;
  char *sequence;
  char *structure;
  int len;
} options;

options move_opt;
  

static void parse_infile(FILE *fp);
int get_list(struct_en*, struct_en*);
int construct_moves_new(const char*, const short*, int , move_str **);
void random_move_pt(short int*);
inline int try_insert_seq(const char*, int, int);
inline int compat(const char, const char);
void mtw_dump_pt(const short*);
gsl_histogram *ini_histogram(const int n,const double min,const double max);

struct_en *list = NULL;
int list_length = 0;
gsl_histogram *h = NULL;

int main() {
  int i,j,count,e,ept;
  short int *pt = NULL;
  short int *ptbak = NULL;
  short int *s0 = NULL;
  short int *s1 = NULL;
  //move_str *mvs = NULL;
  // //char seq[22] = "CCAGUUCCUCACAGGCACAC";
  ////char str[22] = "...(..(((...)))..)..";
  move_opt.INFILE = stdin;
  parse_infile(move_opt.INFILE);
  ////  char *seq = NULL;
  char *str = NULL;

  {
    int hmin=-30.;
    int hmax=50.;
    h = ini_histogram(100,hmin,hmax);
  }
  srand(time(NULL));
  
  pt = vrna_pt_get(move_opt.structure);
  s0 = encode_sequence(move_opt.sequence,0);
  s1 = encode_sequence(move_opt.sequence,1);
  		 
  //mtw_dump_pt(pt);
  {
    e = energy_of_structure(move_opt.sequence,move_opt.structure,0);
    ept = energy_of_structure_pt(move_opt.sequence,pt,s0,s1,0);
    
  }
  //str = vrna_pt_to_db(pt);
  //printf ("structure:\n%s (e=%i) (ept=%i)\n",str,e,ept);
  printf("%s\n", move_opt.sequence);
  print_str(stdout,pt);printf("\n");
  //print_stren(stdout,pt);

  random_move_pt(pt);
  
  gsl_histogram_free(h);
  free(list);
  free(pt);
  free(s0);
  free(s1);
  free(move_opt.sequence);
  free(move_opt.structure);
  return 0;
}

void
random_move_pt(short int *pt)
{
  int i,k,count;
  move_str *mvs=NULL;
  k=0;
  count = construct_moves_new((const char *)move_opt.sequence,pt,1,&mvs);
  //for (i = 0; i<count; i++) {  
  //  printf("%d %d\n", mvs[i].left, mvs[i].right);
  //}
  //printf ("++ applying move #%i: left %i right %i\n",k,mvs[k].left,mvs[k].right);
  
  if(mvs[k].left < 0){
    pt[(int)(fabs(mvs[k].left))] = 0;
    pt[(int)(fabs(mvs[k].right))] = 0;
  }
  else {
    pt[mvs[k].left] = mvs[k].right;
    pt[mvs[k].right] = mvs[k].left;
  }
  print_str(stdout,pt);printf("\n");
  free(mvs);  
}

/**/
static void
parse_infile(FILE *fp)
{
  char *line=NULL;
  
  line = get_line(fp);
  /* skip comment lines */
  while ((*line == '*')||(*line == '\0')||(*line == '>')) {
    free(line);
    line = get_line(fp);
  }
  move_opt.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  move_opt.structure = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(move_opt.sequence != NULL); assert(move_opt.structure != NULL);
  sscanf(line, "%s", move_opt.sequence);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", move_opt.structure);
  free (line);
  move_opt.len = strlen(move_opt.sequence);

}

int
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
          if (try_insert_seq(seq,i,j)) {
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
try_insert_seq(const char *seq,
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

gsl_histogram *ini_histogram(const int n,const double min,const double max) {
  gsl_histogram *a = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(a,min,max);
  return a;
}
