/*
   wl_rna.c RNA-related routines for Wang-Landau sampling
   Last changed Time-stamp: <2014-06-26 15:48:09 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "config.h"
#include "wl_options.h"
#include "globals.h"
#include "wl_rna.h"

static void subopt_first_bin_RNA(float);

/* ==== */
void
initialize_RNA (const char *seq)
{
  model_detailsT md;
  vrna_fold_compound *vc = NULL;
  set_model_details(&md); /* use current global model */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
  s0 = vrna_seq_encode_simple(wanglandau_opt.sequence,&(P->model_details));
  s1 = vrna_seq_encode(wanglandau_opt.sequence,&(P->model_details));

  /* compute mfe */
  mfe = vrna_fold(vc,NULL);
  vrna_free_fold_compound(vc);
  printf ("[[initialize_RNA()]]\nmfe = %6.2f\n",mfe);
}

/* ==== */
void
pre_process_RNA(void)
{
  subopt_first_bin_RNA(wanglandau_opt.erange); 
}

/* ==== */
void
post_process_RNA(void)
{
  return;
}

/* ==== */
static void
subopt_first_bin_RNA(float e)
{
  int strucs, have_lowest_bin;
  size_t i;
  SOLUTION *sol=NULL;

  have_lowest_bin = 0;
  printf("[[subopt_first_bin()]]\ne=%g\n",e);
  sol = subopt(wanglandau_opt.sequence, NULL, e*100, NULL);
  for (strucs = 0; sol[strucs].structure != NULL; strucs++){
    gsl_histogram_find(h,sol[strucs].energy, &i);
    if (i == 0){have_lowest_bin=1;}
    printf("%s %6.2f %d\n",sol[strucs].structure,sol[strucs].energy,i);
    gsl_histogram_increment(h,sol[strucs].energy);
     // TODO: so_structs[bin] += 1; ueber eigenes s histogram
    free(sol[strucs].structure);
  }
  free(sol);
  if (have_lowest_bin != 1){
    fprintf(stderr,"Lowest bin has not been populated in subopt_first_bin()\n");
    fprintf(stderr,"Please decrease the number of bins for this simulation\n");
    fprintf(stderr,"exiting ...\n");
    exit(EXIT_FAILURE);
  }
    
}


