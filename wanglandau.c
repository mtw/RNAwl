/*
  wanglandau.c : main computation routines for Wang-Landau sampling
  Last changed Time-stamp: <2014-06-26 15:49:07 mtw>

  Literature:
  Landau, PD and Tsai, S-H and Exler, M (2004) Am. J. Phys. 72:(10) 1294-1302
  A new approach to Monte Carlo simulations in statistical physics:
  Wang-Landau sampling
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "globals.h"
#include "wl_options.h"
#include "wl_rna.h"
#include "moves.h"

/* functions */
static void initialize_wl(void);
static gsl_histogram *ini_histogram(const int,const double,const double);
static void wl_montecarlo(const char*, char *);

/* variables */
static int iterations=0;   /* # of iterations (modifications with f) */
static long int steps=0;   /* # of WL steps */


/* ==== */
void
wanglandau(void)
{
  initialize_wl();     /* set function pointers for current model;
			  allocate histograms */
  pre_process_model(); /* get normalization factor for histogram by
			  populating the first bin */
  printf("[[wanlandau()]]\n");
  gsl_histogram_fprintf(stderr,h,"%6.3g","%6g");
  wl_montecarlo(wanglandau_opt.sequence, wanglandau_opt.structure);
  post_process_model();
  return;
}

/* ==== */
static void
initialize_wl(void)
{
  double lo,hi, hmin, hmax,erange;
  srand(time(NULL));
  printf("[[initialize_wl()]]\n");
  /* assign function pointers */
  initialize_model = initialize_RNA;  /* for RNA */
  pre_process_model  = pre_process_RNA;
  post_process_model = post_process_RNA;

  initialize_model(wanglandau_opt.sequence); /* set energy paramters for
						current model; get mfe */
  hmin=floor(mfe);
  hmax=5*fabs(mfe);
  //printf ("-> mfe is %4.2f | hmin is %4.2f | hmax is %4.2f <-\n", mfe,hmin,hmax);
  h = ini_histogram(wanglandau_opt.bins,(int)hmin,hmax);
  gsl_histogram_set_ranges_uniform(h,hmin,hmax);
   
  // populate lowest bin with true DOS from subopt
  gsl_histogram_get_range(h,0,&lo,&hi);
  wanglandau_opt.erange=(float)fabs(mfe-hi+0.01);
 
  printf("bin 1: %6.3g -- %6.3g wl_opt.erange=%6.3f\n",lo,hi,wanglandau_opt.erange);
  //gsl_histogram_fprintf(stderr,h,"%6.3g","%6g");
}

/* ==== */
gsl_histogram *
ini_histogram(const int n,
	      const double min,
	      const double max)
{
  gsl_histogram *a = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(a,min,max);
  return a;
}

/* ==== */
static void
wl_montecarlo(const char *seq, char *struc)
{
  short *pt=NULL;
  int e,enew,emove,verbose;
  move_str m;
  verbose = 0;
  
  pt = vrna_pt_get(struc);
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);
  e = vrna_eval_structure_pt(wanglandau_opt.sequence,pt,P);
  printf("%s\n", wanglandau_opt.sequence);
  print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);

  m = get_random_move_pt(seq,pt,verbose);
  emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
  apply_move_pt(pt,m);
  //mtw_dump_pt(pt);
  enew = e + emove;
  //printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
  print_str(stdout,pt);printf(" %6.2f\n",(float)enew/100);
  //e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
  //print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
}

void
sighandler (int signum)
{
  fprintf(stderr, "steps=%12li\n", steps);
    fflush(stderr);
  signal(SIGUSR1, sighandler);
}

void
wanglandau_free_memory(void)
{
  gsl_histogram_free(h);
  free(pt);
  free(s0);
  free(s1);
  free(wanglandau_opt.sequence);
  free(wanglandau_opt.structure);
}
