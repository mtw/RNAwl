#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

#define MAXG 2048

static gsl_histogram *ini_histogram(const int,const double,const double);
static short histogram_is_flat(const gsl_histogram *);

int bins = 10;

int main (void)
{
  gsl_histogram *m = NULL;
  short is_flat = 0;
  m = ini_histogram(bins,0,100);
  gsl_histogram_increment(m,5);
  gsl_histogram_increment(m,25);
  gsl_histogram_increment(m,55);
  gsl_histogram_increment(m,57);
  gsl_histogram_increment(m,65);
  gsl_histogram_fprintf(stderr,m,"%6.3g","%6g");

  if (histogram_is_flat(m)){
    printf("Histogram is flat\n");
  }
  else {printf ("Histogram is NOT flat\n");}
  return;
}

/* ==== */
static short
histogram_is_flat(const gsl_histogram *z)
{
  double val,avg,sum = 0.0,flat = 0.8;
  size_t lbin,gbin;        /* lowest/greatest populated bin */
  int i,b=0,is_flat=1;

  fprintf(stderr,"[[histogram_is_flat()]]\n");
  for (i=0; i<bins; i++){
    val = gsl_histogram_get(z,i);
    if (val>0){
      lbin=i;
      break;
    }
  }
  
  printf("lbin is %i\n",lbin);
  for (i=bins-1; i>=0; i--){
    val = gsl_histogram_get(z,i);
    if (val > 0){
      gbin=i;
      break;
    }
  }
  printf("gbin is %i\n",gbin);

  // compute average over histogram
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
      if( val < flat*avg){
	is_flat = 0;
	break;
      }
    }
  }

  
 
  return is_flat;
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
