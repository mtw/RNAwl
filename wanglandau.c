/*
  wanglandau.c : main computation routines for Wang-Landau sampling
  Last changed Time-stamp: <2014-07-31 22:42:29 mtw>

  Literature:
  Landau, PD and Tsai, S-H and Exler, M (2004) Am. J. Phys. 72:(10) 1294-1302
  A new approach to Monte Carlo simulations in statistical physics:
  Wang-Landau sampling
*/

/*
  TODO
  - ensure trudosbins are not scales before output
  - check ln(g) in truedosbins -> must stay constant
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include <sys/time.h>
#include <assert.h>
#include "globals.h"
#include "wl_options.h"
#include "wl_rna.h"
#include "moves.h"
#include <gsl/gsl_rng.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#else
#include <time.h>
#endif

#define MIN2(A, B)  ((A) < (B) ? (A) : (B))
#define MAX2(A, B)  ((A) > (B) ? (A) : (B))

/* functions */
static void initialize_wl(void);
static void initialize_dos_estimate(void);
static void wl_montecarlo(char *);
static gsl_histogram * scale_dos(gsl_histogram *);
static void output_dos(const gsl_histogram *, const char);
static short histogram_is_flat(const gsl_histogram *);
static double partition_function(const gsl_histogram *);
static gsl_histogram *ini_histogram_uniform(const int,const double,const double);

/* variables */
static int iterations = 0;    /* #iterations (modifications with f) */
static int maxbin = -1;       /* index of highest bin */
static unsigned long steps = 0; /* # of WL steps */
static unsigned long seed;    /* random seed */
static gsl_rng *r = NULL;     /* GSL random number generator */
static struct timespec ts;    /* timespec struct for random seed */
static double rnum;           /* random number */

/* arrays */
static gsl_histogram *g = NULL;  /* DoS histogram */
static char *out_prefix=NULL;    /* prefix for output */

/* ==== */
void
wanglandau(void)
{
  initialize_wl();           /* set function pointers for current
				model; allocate histograms */
  pre_process_model();       /* get normalization factor for histogram
				by populating the first bin */
  initialize_dos_estimate();  /* set initial DOS estimate to start
				 with */
  wl_montecarlo(wanglandau_opt.structure);
  // scale_normalize_DOS();
  post_process_model();
  return;
}

/* ==== */
static void
initialize_wl(void)
{
  int bins,fnlen=512;
  char *res_string=NULL;
  double low,high,lo,hi,hmin,hmax,erange,*range = NULL;
  
  srand(time(NULL));
  if(wanglandau_opt.verbose){
    printf("[[initialize_wl()]]\n");
  }
  /* assign function pointers */
  initialize_model = initialize_RNA;  /* for RNA */
  pre_process_model  = pre_process_RNA;
  post_process_model = post_process_RNA;
  
  /* set energy paramters for current model; compute mfe */
  initialize_model(wanglandau_opt.sequence); 

  range = (double*)calloc((wanglandau_opt.bins+1), sizeof(double));
  assert(range!=NULL);

  /* initialize histograms */
  hmin=mfe;
  if(wanglandau_opt.res_given){ /* determine histogram ranges manually */
    int i;
    range[0]=mfe;
    for(i=1;i<=wanglandau_opt.bins;i++){
      range[i]=range[i-1]+wanglandau_opt.res;
    }
    if(wanglandau_opt.verbose){
      /* info output */
      fprintf(stderr,"#allocating %lu bins of width %g\n",
	      wanglandau_opt.bins,wanglandau_opt.res);
      /* fprintf(stderr,"#histogram ranges:\n #");
	 for(i=0;i<=wanglandau_opt.bins;i++){
	 fprintf(stderr, "%6.2f ",range[i]);
      
      }
      fprintf(stderr,"\n");
      */
    }
  
    if(wanglandau_opt.max_given){
      hmax=wanglandau_opt.max;
    }
    else{
      wanglandau_opt.max=range[wanglandau_opt.bins]; /* the last element */
      hmax=wanglandau_opt.max;
    }
    h = gsl_histogram_alloc(wanglandau_opt.bins);
    g = gsl_histogram_alloc(wanglandau_opt.bins);
    s = gsl_histogram_alloc(wanglandau_opt.bins);
    gsl_histogram_set_ranges(h,range,(wanglandau_opt.bins+1));
    gsl_histogram_set_ranges(g,range,(wanglandau_opt.bins+1));
    gsl_histogram_set_ranges(s,range,(wanglandau_opt.bins+1));
  }
  else{  /* determine histogram ranges automatically */
    if(wanglandau_opt.max_given){
      hmax = wanglandau_opt.max;
    }
    else{
      hmax=20*fabs(mfe);
    }
    h = ini_histogram_uniform(wanglandau_opt.bins,hmin,hmax);
    g = ini_histogram_uniform(wanglandau_opt.bins,hmin,hmax);
    s = ini_histogram_uniform(wanglandau_opt.bins,hmin,hmax);
  }
  fprintf (stderr, "# sampling energy range is %6.2f - %6.2f\n",
	   hmin,hmax);
  free(range);
  
  /* get the energy range up to which we will compute true DOS via
     RNAsubopt */
  gsl_histogram_get_range(s,0,&lo,&hi);
  low = lo;
  gsl_histogram_get_range(s,(wanglandau_opt.truedosbins-1),&lo,&hi);
  high = hi;
  wanglandau_opt.erange=(float)fabs(mfe-high+0.01);
  if(wanglandau_opt.verbose){
    printf("Using true DOS for bins 0-%d: (%6.3g -- %6.3g) wl_opt.erange=%6.3f\n",
	   (wanglandau_opt.truedosbins-1),low,high,wanglandau_opt.erange);
  }
  
  /* prepare gsl random-number generation */
  (void) clock_gettime(CLOCK_REALTIME, &ts);
  if(wanglandau_opt.seed_given){
    seed = wanglandau_opt.seed;
  }
  else {
    seed =   ts.tv_sec ^ ts.tv_nsec;
  }
  fprintf(stderr, "initializing random seed: %d\n",seed);
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set( r, seed );
  /* end gsl */

  /* make prefix for output */
  out_prefix = (char*)calloc(fnlen, sizeof(char));
  res_string = (char*)calloc(16, sizeof(char));
  sprintf(res_string,"%3.1f", wanglandau_opt.res);
  strcpy(out_prefix, wanglandau_opt.basename); strcat(out_prefix, ".res");
  strcat(out_prefix, res_string); strcat(out_prefix, ".");
  free(res_string);
  return;
}

/* ==== */
static void
initialize_dos_estimate(void)
{
  int i;
  const size_t n = s->n; /* nr of bins */

  if(wanglandau_opt.verbose){
    fprintf(stderr, "[[initialize_dos_estimate()]]\n");
  }
  
  if (wanglandau_opt.truedosbins_given){
    /* initialize the first n bins of g with true DOS as computed by
       RNAsubopt */
    if(wanglandau_opt.verbose){
      fprintf(stderr, "initializing the lowest %d bins of DOS estimate g with true DOS values from subopt:\n",
	      wanglandau_opt.truedosbins);
    }
    for (i=0;i<wanglandau_opt.truedosbins;i++){
      g->bin[i]=log(s->bin[i]);  /* get corresponding value from s histogram */ }
    for(i=wanglandau_opt.truedosbins;i<n;i++){
      g->bin[i]=g->bin[wanglandau_opt.truedosbins-1];
    }
  }
  else{
    /* initialize g uniformly with 0 */
    if (wanglandau_opt.verbose){
      fprintf(stderr, "initializing DOS estimate g with zeros:\n");
    }
  }
  if (wanglandau_opt.verbose){
    gsl_histogram_fprintf(stderr,g,"%6.2f","%30.6f");
    fprintf(stderr,"+++\ndone initializing\n");
  }
}



/* ==== */
static void
wl_montecarlo(char *struc)
{
  short *pt=NULL;
  move_str m;
  int e,enew,emove,eval_me,status,debug=1;
  long int crosscheck=1000000; /* used for convergence checks */
  long int crosscheck_limit = 100000000000000000;
  double g_b1,g_b2,prob,lnf = 1.;  /* log modification parameter f */
  size_t b1,b2;                    /* indices in g/h corresponding to
				      old/new energies */
  gsl_histogram *gcp=NULL; /* clone of g used during crosscheck output */ 

  eval_me = 1; /* paranoid checking of neighbors against RNAeval */
  if (wanglandau_opt.verbose){
    printf("[[wl_montecarlo()]]\n");
  }
  pt = vrna_ptable(struc);
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.temperature = wanglandau_opt.T;
  vrna_fold_compound_t *vc = vrna_fold_compound(wanglandau_opt.sequence,&md,VRNA_OPTION_EVAL_ONLY);
  e = vrna_eval_structure_pt(vc,pt);
  
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
  print_str(stderr,pt);
  printf(" (%6.2f) bin:%d\n",(float)e/100,b1);
  if (wanglandau_opt.verbose){
    fprintf(stderr,"\nStarting MC loop ...\n");
  }
  while (lnf > wanglandau_opt.ffinal) {
    if(wanglandau_opt.debug){
      fprintf(stderr,"\n==================\n");
      fprintf(stderr,"in while: lnf=%8.6f\n",lnf);
      fprintf(stderr,"steps: %d\n",steps);
      fprintf(stderr,"current histogram g:\n");
      gsl_histogram_fprintf(stderr,g,"%6.2f","%30.6f");
      fprintf(stderr,"\n");
      print_str(stderr,pt);
      fprintf(stderr, " (%6.2f) bin:%d\n",(float)e/100,b1);
      /*  mtw_dump_pt(pt); */
    }
    /* make a random move */
    m = get_random_move_pt(wanglandau_opt.sequence,pt);
    /* compute energy difference for this move */
    emove = vrna_eval_move_pt(vc,pt,m.left,m.right);
    /* evaluate energy of the new structure */
    enew = e + emove;
    if(wanglandau_opt.debug){
      fprintf(stderr,
	      "random move: left %i right %i enew(%6.4f)=e(%6.4f)+emove(%6.4f)\n",
	      m.left,m.right,(float)enew/100,(float)e/100,(float)emove/100);
    }

    /* ensure the new energy is within sampling range */
    if ((float)enew/100 >= wanglandau_opt.max){
      fprintf(stderr,
	      "New structure has energy %6.2f >= %6.2f (upper energy bound)\n",
	      (float)enew/100,wanglandau_opt.max);
      fprintf(stderr,"Please increase --bins or adjust --max! Exiting ...\n");
      exit(EXIT_FAILURE);
    }
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
    g_b1 = gsl_histogram_get(g,b1);
    g_b2 = gsl_histogram_get(g,b2);
      
    /* core MC steps */
    prob = MIN2(exp(g_b1 - g_b2), 1.0);
    rnum =  gsl_rng_uniform (r);
    
    if ((prob == 1 || (rnum <= prob)) ) { /* accept & apply the move */
      apply_move_pt(pt,m);
      if(wanglandau_opt.debug){
	print_str(stderr,pt);
	fprintf(stderr, " %6.2f bin:%d [A]\n", (float)enew/100,b2);
      }
      b1 = b2;
      e = enew;
    }
    else { /* reject the move */
      if(wanglandau_opt.debug){
	print_str(stderr,pt);
	fprintf(stderr, " (%6.2f) bin:%d [R]\n", (float)enew/100,b2);
       }
    }
    
    /* update histograms g and h */
    if(wanglandau_opt.truedosbins_given && b2 <= wanglandau_opt.truedosbins){
      /* do not update if b2 <= truedosbins, i.e. keep true DOS values
	 in those bins */
      if (wanglandau_opt.debug){
	fprintf(stderr, "NOT UPDATING bin %d\n",b1);
      }
    } else{
      if(wanglandau_opt.debug){
	fprintf(stderr, "UPDATING bin %d\n",b1); 
      }
      status = gsl_histogram_increment(h,(float)e/100);
      status = gsl_histogram_accumulate(g,(float)e/100,lnf);
    }
    maxbin = MAX2(maxbin,(int)b1);
   
    // stuff that can be skipped 
    /*
      printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
      print_str(stderr,pt);printf(" %6.2f bin:%d\n",(float)enew/100,b2);
      e = vrna_eval_structure_pt(wanglandau_opt.sequence,pt,P);
      if (eval_me == 1 && e != enew){
      fprintf(stderr, "energy evaluation against vrna_eval_structure_pt() mismatch... HAVE %6.2f != %6.2f (SHOULD BE)\n",(float)enew/100, (float)e/100);
      exit(EXIT_FAILURE);
      }
      print_str(stderr,pt);printf(" %6.2f\n",(float)e/100);
    */
    // end of stuff that can be skipped

    /* output DoS every x*10^(1/4) steps, starting with x=10^6 (we
       used this fopr comparing perfomance and convergence of
       different DoS sampling methods */
    if((steps % crosscheck == 0) && (crosscheck <= crosscheck_limit)){
      fprintf(stderr,"# crosscheck reached %li steps ",crosscheck);
      gcp = gsl_histogram_clone(g);
      if(wanglandau_opt.verbose){
	fprintf(stderr,"## gcp before scaling\n");
	gsl_histogram_fprintf(stderr,gcp,"%6.2f","%30.6f");
      }
      scale_dos(gcp); /* scale estimated g; make ln(g[0])=0 */
      if(wanglandau_opt.verbose){
	fprintf(stderr,"## gcp after scaling\n");
	gsl_histogram_fprintf(stderr,gcp,"%6.2f","%30.6f");
      }
      double Z = partition_function(gcp);
      output_dos(gcp,'s');
      fprintf(stderr, "Z=%10.4g\n", Z);
      crosscheck *= (pow(10, 1.0/4.0));
      gsl_histogram_free(gcp);
      fprintf(stderr,"->  new crosscheck will be performed at %li steps\n", crosscheck);
    }
    
    if(steps % wanglandau_opt.checksteps == 0) {
      if( histogram_is_flat(h) ) {
	lnf /= 2;
	fprintf(stderr,"# steps=%20li | f=%12g | histogram is FLAT\n",
		steps,lnf);
	gsl_histogram_reset(h);
      }
      else {
	fprintf(stderr, "# steps=%20li | f=%12g | histogram is NOT FLAT\n",
		steps,lnf);
      }
      output_dos(g,'l');
    }
    
    /* stop criterion */
    if(steps % wanglandau_opt.steplimit == 0){
      fprintf(stderr,"maximun number of MC steps (%li) reached, exiting ...",
	      wanglandau_opt.steplimit);
      break;
    }

  } /* end while */

  vrna_fold_compound_free(vc);
  free(pt); 
  return;
}


/* ==== */
static gsl_histogram *
ini_histogram_uniform(const int n,
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
  size_t lbin,gbin;         /* lowest/greatest populated bin */
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
/* Z = \sum{E} g(e)*e^{-E-kT} */
static double
partition_function(const gsl_histogram *y)
{
  int i;
  float T;
  double kT,Z=0.;
  const size_t n = y->n; /* nr of bins */

  T = 273.15 + wanglandau_opt.T;
  kT  = 0.00198717*4.16*T;
  
  for(i=0;i<n;i++){
    Z += y->bin[i] * exp(-1*y->bin[i]/kT);
  }
  return Z;
}

/* ==== */
static gsl_histogram *
scale_dos(gsl_histogram *y)
{
  int i,maxbin;
  size_t bins;
  double  maxval=-1., sum=0., x=0, factor=0., GZero=0,  exp_G_norm=0.;
  const size_t n = y->n; /* nr of bins */
  
  /* FIRST: scale it via the ground state */
  /* ln[gn(E)] = ln[g(E)]-ln[g(Egs)]+ln[Q] */
  /* where Q is the # of structures found in the lowest bin/groundstate */

  bins = gsl_histogram_bins(y);
  /* get value of y[0] */
  GZero = gsl_histogram_get(y,0);
  /* compute scaling factor just from the very lowest bin for now */
  for(i=0;i<1;i++){
    factor += gsl_histogram_get(s,i);
  }
  /* subtract g[0] [ln(g(Egs))] from each entry to get smaller numbers
     and add scaling factor*/
  for (i=wanglandau_opt.truedosbins-1;i<n;i++){
    if(y->bin[i] == 0.){ continue;}
    else{y->bin[i] += log(factor)-GZero;}
  }
  
  /* exponentiate to get effective DOS */
  /*
  for(i=0;i<n;i++){
    y->bin[i]=exp(y->bin[i]);
  }
  */
  
  output_dos(y,'s');  
  
  return y;
}

/* ==== */
static void
output_dos(const gsl_histogram *x, const char T)
{
  int i,fnlen;
  FILE *dos_fp=NULL;
  char *dos_fn=NULL, *lDoS_suffix="lDoS", *sDoS_suffix="sDoS";
  char s[50];
  double val,lo,hi;
 
  
  sprintf(s,"%li",steps);
  fnlen = strlen(out_prefix)+strlen(lDoS_suffix)+64;
  dos_fn = (char*)calloc(fnlen,sizeof(char));
  strcpy(dos_fn, out_prefix);
  strcat(dos_fn, s);
  strcat(dos_fn,".");
  
  switch (T){
  case 'l':  /* output logarithmic g, usually during the calculation  */
    strcat(dos_fn, lDoS_suffix);
    break;
  case 's':  /* output scaled g, eg for in-process convergence checks */
    strcat(dos_fn, sDoS_suffix);
    break;
  default:
    fprintf (stderr, "%s:%d output_dos(): No handler for type %c",
	     __FILE__, __LINE__, T);
    exit(EXIT_FAILURE);
  }
  
  dos_fp = fopen(dos_fn, "w+");
  fprintf(dos_fp, "# estimated DOS after %li steps\n",steps);
  fprintf(dos_fp, "# sampling range: %6.2f -- %6.2f\n",
	  gsl_histogram_min(g),gsl_histogram_max(g));
  fprintf(dos_fp, "# bin resolution: %g\n",wanglandau_opt.res);

  /* loop over histogram g */
  for (i=0;i<=maxbin;i++){
    val = gsl_histogram_get(x,i);
    if (val == 0.){continue;}
    gsl_histogram_get_range(x,i,&lo,&hi);
    fprintf(dos_fp,"%6.2f\t%20.6f\n",lo+(hi-lo)/2,val);
  }  
  fclose(dos_fp);
  free(dos_fn);
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
  free(wanglandau_opt.sequence);
  free(wanglandau_opt.structure);
  free(wanglandau_opt.basename);
  free(out_prefix);
  dealloc_gengetopt();
  return;
}
