#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "moves.h"
#include <gsl/gsl_histogram.h>
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/subopt.h>

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
move_str get_random_move_pt(const short int*);
void apply_move_pt(short int *,const move_str);
void subopt_first_bin(double);
inline int try_insert_seq(const char*, int, int);
inline int compat(const char, const char);
void mtw_dump_pt(const short*);
void ini_RNA(const char*);
gsl_histogram *ini_histogram(const int n,const double min,const double max);

struct_en *list = NULL;
int list_length = 0;
gsl_histogram *h = NULL;
model_detailsT md;
paramT *P = NULL;
vrna_fold_compound *vc = NULL;

int main() {
  int e,enew,emove;
  float mfe;
  short int *pt,*s0,*s1;
  move_str m;
  double erange;
  
  move_opt.INFILE = stdin;
  parse_infile(move_opt.INFILE);
  ini_RNA(move_opt.sequence);
  srand(time(NULL));
  
  { // compute mfe
    mfe = vrna_fold(vc,NULL);
    destroy_fold_compound(vc);
    printf ("mfe = %6.2f\n",mfe);
  }
  { // alloc and initalize histogram
    double lo,hi, hmin=floor(mfe);
    int hmax=5*fabs(mfe);
    h = ini_histogram(20,(int)hmin,hmax);
    gsl_histogram_set_ranges_uniform(h,hmin,hmax);
  
    //gsl_histogram_increment(h,mfe);
    gsl_histogram_fprintf(stderr,h,"%6.2g","%6g");
   
  }
  { // populate lowest bin with true DOS from subopt
    gsl_histogram_get_range(h,0,&lo,&hi);
    erange=fabs(mfe-hi+0.01);
    printf("bin 1: %g -- %g; running subopt -e %g\n",lo,hi,erange);
    subopt_first_bin(erange);
    gsl_histogram_fprintf(stderr,h,"%6.2g","%6g");
  }
  
  pt = vrna_pt_get(move_opt.structure);
  s0 = get_sequence_encoding(move_opt.sequence,0,&(P->model_details));
  s1 = get_sequence_encoding(move_opt.sequence,1,&(P->model_details));
  		 
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);
  e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
  printf("%s\n", move_opt.sequence);
  print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);

  m = get_random_move_pt(pt);
  emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
  apply_move_pt(pt,m);
  //mtw_dump_pt(pt);
  enew = e + emove;
  //printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
  print_str(stdout,pt);printf(" %6.2f\n",(float)enew/100);
  //e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
  //print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
  
  gsl_histogram_free(h);
  free(list);
  free(pt);
  free(s0);
  free(s1);
  free(move_opt.sequence);
  free(move_opt.structure);
  return 0;
}

/*
  compute a random move on a pair table
  returns move operations to be applied to pt in order to perform the move
 */
move_str
get_random_move_pt(const short int *pt)
{
  move_str r,*mvs=NULL;
  
  int count = construct_moves_new((const char *)move_opt.sequence,pt,1,&mvs);

  /*
    for (i = 0; i<count; i++) {  
    printf("%d %d\n", mvs[i].left, mvs[i].right);
    }
    printf ("++ applying move #%i: left %i right %i\n",i,mvs[i].left,mvs[i].right);
  */
  
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

void
ini_RNA (const char *seq)
{
  set_model_details(&md); /* use current global model */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
}

void
subopt_first_bin(double e)
{
  int strucs;
  SOLUTION *sol=NULL;
  printf("in subopt_first_bin: e=%g\n",e);
  sol = subopt(move_opt.sequence, NULL, e*100, NULL);
  for (strucs = 0; sol[strucs].structure != NULL; strucs++){
    printf("%s %6.2f\n",sol[strucs].structure,sol[strucs].energy);
    gsl_histogram_increment(h,sol[strucs].energy);
    // TODO: so_structs[bin] += 1; ueber eigenes s histogram
    free(sol[strucs].structure);
  }
  free(sol);
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
