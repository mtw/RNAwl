/*
  wanglandau.c : main computation routines for Wang-Landau sampling
  Last changed Time-stamp: <2014-07-03 00:27:55 mtw>

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
#include <sys/time.h>
#include <assert.h>
#include "globals.h"
#include "wl_options.h"
#include "wl_rna.h"
#include "moves.h"
#include <gsl/gsl_rng.h>

#define MIN2(A, B)  ((A) < (B) ? (A) : (B))
#define MAX2(A, B)  ((A) > (B) ? (A) : (B))

/* functions */
static void initialize_wl(void);
static void wl_montecarlo(char *);
static void scale_normalize_DOS(void);
static void output_dos(double *dos, char T);
static short histogram_is_flat(const gsl_histogram *);
static gsl_histogram *ini_histogram(const int,const double,const double);

/* variables */
static int iterations = 0;       /* # of iterations (modifications with f) */
static int maxbin = -1;          /* index of highest bin */
static unsigned long steps = 0;  /* # of WL steps */
static unsigned long seed;       /* random seed */
static gsl_rng *r = NULL;        /* GSL random number generator */
static struct timespec ts;       /* timespec struct needed for random seed */
static double rnum;              /* random number */

/* arrays */
gsl_histogram *g = NULL;         /* DoS histogram */
static char *out_prefix=NULL;    /* prefix for output */

/* ==== */
void
wanglandau(void)
{
  initialize_wl();     /* set function pointers for current model;
			  allocate histograms */
  pre_process_model(); /* get normalization factor for histogram by
			  populating the first bin */
  
 
  if(wanglandau_opt.verbose){
    printf("[[wanlandau()]]\n");
  }
  wl_montecarlo(wanglandau_opt.structure);
  scale_normalize_DOS();
  post_process_model();
  return;
}

/* ==== */
static void
initialize_wl(void)
{
  int bins;
  double low,high,lo,hi,hmin,hmax,erange;
  
  srand(time(NULL));
  if(wanglandau_opt.verbose){
    printf("[[initialize_wl()]]\n");
  }
  /* assign function pointers */
  initialize_model = initialize_RNA;  /* for RNA */
  pre_process_model  = pre_process_RNA;
  post_process_model = post_process_RNA;

  initialize_model(wanglandau_opt.sequence); /* set energy paramters for
						current model; get mfe */
  //hmin=floor(mfe);
  hmin=mfe-0.10001;
  hmax=5*fabs(mfe);
  /*printf ("-> mfe is %4.2f | hmin is %4.2f | hmax is %4.2f <-\n",
    mfe,hmin,hmax); */

  /* prepare histogram h (holds energy levels seen in the current
     iteration) */
  h = ini_histogram(wanglandau_opt.bins,(int)hmin,hmax);
  bins = gsl_histogram_bins(h);
  if(wanglandau_opt.verbose){
    fprintf(stderr, "histogram h allocated with %zu bins\n",bins);
  }
 
  /* prepare histogram g (holds estimation for DOS */
  g = ini_histogram(wanglandau_opt.bins,(int)hmin,hmax);
  bins = gsl_histogram_bins(g);
  if(wanglandau_opt.verbose){
    fprintf(stderr, "histogram g allocated with %d bins\n",bins);
  }

  /* prepare histogram s (required for normalization) */
  s = ini_histogram(wanglandau_opt.bins,(int)hmin,hmax);
  bins = gsl_histogram_bins(s);
  if(wanglandau_opt.verbose){
    fprintf(stderr, "histogram s allocated with %d bins\n",bins);
  }

  dos = (double*)calloc(wanglandau_opt.bins,sizeof(double));
  assert(dos != NULL);
  

  /* get the energy range up to which we will compute true DOS via
     RNAsubopt, which will be required later for normalization */
  gsl_histogram_get_range(s,0,&lo,&hi);
  low = lo;
  gsl_histogram_get_range(s,(wanglandau_opt.norm-1),&lo,&hi);
  high = hi;
  wanglandau_opt.erange=(float)fabs(mfe-high+0.01);
  if(wanglandau_opt.verbose){
    printf("bins 0-%d: %6.3g -- %6.3g wl_opt.erange=%6.3f\n",
	   (wanglandau_opt.norm-1),low,high,wanglandau_opt.erange);
    gsl_histogram_fprintf(stderr,h,"%6.3g","%6g");
  }
  
  /* prepare gsl random-number generation */
  (void) clock_gettime(CLOCK_REALTIME, &ts);
  if(wanglandau_opt.seed_given){
    seed = wanglandau_opt.seed;
  }
  else {
    seed =   ts.tv_sec ^ ts.tv_nsec;
  }
  printf("initializing random seed: %d\n",seed);
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set( r, seed );

  return;
}

/* ==== */
static void
wl_montecarlo(char *struc)
{
  short *pt=NULL;
  int e,enew,emove,eval_me,status,debug=1;
  double g_b1,g_b2,prob,lnf = 1;   /* logarithmic modification parameter f */
  size_t b1,b2;                    /* indices in g/h corresponding to old/new energies */
  move_str m;

  eval_me = 1; /* paranoid checking of neighbors against RNAeval */
  if (wanglandau_opt.verbose){
    printf("[[wl_montecarlo()]]\n");
  }
  pt = vrna_pt_get(struc);
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);
  e = vrna_eval_structure_pt(wanglandau_opt.sequence,pt,P);
  
  /* determine bin where the start structure goes */
  status = gsl_histogram_find(g,(float)e/100,&b1);
  if (status) {
    if (status == GSL_EDOM){
      printf ("error: %s\n", gsl_strerror (status));
    }
    else {fprintf(stderr, "GSL error: gsl_errno=%d\n",status);}
    exit(EXIT_FAILURE);
  }
  printf("%s\n", wanglandau_opt.sequence);
  print_str(stdout,pt);printf(" %6.2f bin:%d\n",(float)e/100,b1);
  if (wanglandau_opt.verbose){
    fprintf(stderr,"Starting MC loop ...\n");
  }
  while (lnf > wanglandau_opt.ffinal) {
    /* make a random move */
    m = get_random_move_pt(wanglandau_opt.sequence,pt,wanglandau_opt.verbose);
    /* compute energy difference for this move */
    emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
    /* evaluate energy of the new structure */
    enew = e + emove;
    /* determine bin where the new structure goes */
    status = gsl_histogram_find(g,(float)enew/100,&b2);
    if (status) {
      if (status == GSL_EDOM){
	printf ("error: %s\n", gsl_strerror (status));
      }
      else {fprintf(stderr, "GSL error: gsl_errno=%d\n",status);}
      exit(EXIT_FAILURE);
    }
    steps++;  /* # of MC steps performed so far */
    /* lookup current values for bins b1 and b2 */
    if (debug == 1){
      fprintf(stderr,"==================\n");
      gsl_histogram_fprintf(stderr,g,"%6.3g","%8.6g");
      fprintf(stderr,"\n");
    }
    g_b1 = gsl_histogram_get(g,b1);
    g_b2 = gsl_histogram_get(g,b2);
    if(debug ==1){
      fprintf(stderr,"b1=%i g_b1: %6.2f | b2=%i g_b2: %6.2f\n",b1,g_b1,b2,g_b2);
    }
    
    /* core MC steps */
    prob = MIN2(exp(g_b1 - g_b2), 1.0);
    rnum =  gsl_rng_uniform (r);
    /* if(wanglandau_opt.verbose){ printf ("rnum is %.5f\n", rnum); }*/
    
    if ((prob == 1 || (rnum <= prob)) ) { /* accept the move */
      /* apply the move */
      apply_move_pt(pt,m);
      /* fprintf(stderr, "accepting %s %6.4f\n", neighbor, e); */
      //mtw_dump_pt(pt);
      b1 = b2;
      e = enew;
    }
    else { /* reject the move */
      ;
    }
    /* update histogram h */
    status = gsl_histogram_increment(h,(float)e/100);
    if (status) {
      if (status == GSL_EDOM){
	printf ("increment h error: %s\n", gsl_strerror (status));
	printf ("while trying to increment bin at value %d\n",e);
      }
      else {fprintf(stderr, "GSL error: gsl_errno=%d\n",status);}
      exit(EXIT_FAILURE);
    }
    /* update histogram g */
    status = gsl_histogram_accumulate(g,(float)e/100,lnf);
    if (status) {
      if (status == GSL_EDOM){
	printf ("increment g error: %s\n", gsl_strerror (status));
      }
      else {fprintf(stderr, "GSL error: gsl_errno=%d\n",status);}
      exit(EXIT_FAILURE);
    }
    maxbin = MAX2(maxbin,b1);
    
    // stuff that can be skipped 
    /*
      printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
      print_str(stdout,pt);printf(" %6.2f bin:%d\n",(float)enew/100,b2);
      e = vrna_eval_structure_pt(wanglandau_opt.sequence,pt,P);
      if (eval_me == 1 && e != enew){
      fprintf(stderr, "energy evaluation against vrna_eval_structure_pt() mismatch... HAVE %6.2f != %6.2f (SHOULD BE)\n",(float)enew/100, (float)e/100);
      exit(EXIT_FAILURE);
      }
      print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
    */
    // end of stuff that can be skipped
    
    if(steps % (wanglandau_opt.steps) == 0) {
      if( histogram_is_flat(h) ) {
	lnf /= 2;
	fprintf(stderr,"# histogram is flat  f=%12g steps=%20li\n",lnf,steps);
	gsl_histogram_reset(h);
      }
      else {
	fprintf(stderr, "#lnf=%12g steps=%20li not flat\n", lnf, steps);
      }
    }
  } /* end while */
  free(pt); 
  return;
}


/* ==== */
static gsl_histogram *
ini_histogram(const int n,
	      const double min,
	      const double max)
{
  gsl_histogram *a = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(a,min,max);
  return a;
}

/* ==== */
static short
histogram_is_flat(const gsl_histogram *z)
{
  double val,avg,sum = 0.0;
  size_t lbin,gbin;        /* lowest/greatest populated bin */
  int i,b=0,is_flat=1;

  /* get lowest populated bin */
  for (i=0; i<wanglandau_opt.bins; i++){
    val = gsl_histogram_get(z,i);
    if (val>0){
      lbin=i;
      break;
    }
  }
  /* printf("lbin is %i\n",lbin);*/

  // get highest populated bin
  for (i=wanglandau_opt.bins-1; i>=0; i--){
    val = gsl_histogram_get(z,i);
    if (val > 0){
      gbin=i;
      break;
    }
  }
  /* printf("gbin is %i\n",gbin);*/

  // compute average over interval [lbin;gbin]
  for (i=lbin;i<=gbin;i++){
    if ((val = gsl_histogram_get(z,i)) != 0){
      sum += val;
      b++;
    }
  }
  avg = sum/b;

  // evaluate if populated bins are all within accepted range
  for (i=lbin;i<=gbin;i++){
    if ((val = gsl_histogram_get(z,i)) != 0){
      if( val < wanglandau_opt.flat*avg){
	is_flat = 0;
	break;
      }
    }
  }
  return is_flat;
}

/* ==== */
static void
scale_normalize_DOS(void)
{
  int i,maxbin;
  size_t bins;
  double  maxval=-1., sum=0., x=0, factor=0., GZero=0,  exp_G_norm=0.;

  /* FIRST: scale g */
  /* scale it via the ground state */
  /* ln[gn(E)] = ln[g(E)]-ln[g(Egs)]+ln[Q] */
  /* where Q is the # of structures found in the lowest bin/groundstate */

  bins = gsl_histogram_bins(g);
  /* get value of g[0] */
  GZero = gsl_histogram_get(g,0);
  for(i=0;i<wanglandau_opt.norm;i++){
    double val = gsl_histogram_get(s,i);
    factor += val;
  }
  /* subtract g[0] from each entry to get smaller numbers */
  gsl_histogram_shift(g,(-1*GZero));

  for(i=0;i<wanglandau_opt.norm;i++){
    double val = gsl_histogram_get(g,i);
    exp_G_norm += exp(val);
  }
  for(i=0;i<bins;i++){
    double val = gsl_histogram_get(g,i);
    double tmp = val+log(factor)-log(exp_G_norm);
    dos[i] = tmp;
  }

  /* output dos */
  //output_dos(dos, 'DOS');

  /*
    for(i=0;i<wanglandau_opt.bins;i++){
    fprintf(stderr, %6.2
  }
  */
  gsl_histogram_fprintf(stderr,g,"%6.3g","%6g");
  
  /* SECOND: normalize g */
}

/* ==== */
static void
output_dos(double *dos, char T)
{
  return;
}

/* ==== */
void
sighandler (int signum)
{
  fprintf(stderr, "steps=%12li\n", steps);
    fflush(stderr);
  signal(SIGUSR1, sighandler);
  return;
}

/* ==== */
void
wanglandau_free_memory(void)
{
  gsl_histogram_free(h);
  gsl_histogram_free(g);
  gsl_rng_free(r);
  free(pt);
  free(s0);
  free(s1);
  free(dos);
  free(wanglandau_opt.sequence);
  free(wanglandau_opt.structure);
  free(wanglandau_opt.basename);
  free(P);
  dealloc_gengetopt();
  return;
}
