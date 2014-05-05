/* wl_rna.c RNA-related routines for Wang-Landau sampling */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
//#include "options.h"
#include "globals.h"
#include "wl_rna.h"


/* ViennaRNA-related */
paramT *P;
short int *pt,*s0,*s1;

void
initialize_RNA (const char *seq)
{
  model_detailsT md;
  vrna_fold_compound *vc = NULL;
  set_model_details(&md); /* use current global model */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
  s0 = get_sequence_encoding(wanglandau_opt.sequence,0,&(P->model_details));
  s1 = get_sequence_encoding(wanglandau_opt.sequence,1,&(P->model_details));

  /* compute mfe */
  mfe = vrna_fold(vc,NULL);
  destroy_fold_compound(vc);
  printf ("[[initialize_RNA]]: mfe = %6.2f\n",mfe);
}

void
pre_process_RNA(void)
{
  subopt_first_bin_RNA(wanglandau_opt.erange);
  
}

void
subopt_first_bin_RNA(float e)
{
  int strucs;
  SOLUTION *sol=NULL;
  printf("[[subopt_first_bin]]: e=%g\n",e);
  sol = subopt(move_opt.sequence, NULL, e*100, NULL);
  for (strucs = 0; sol[strucs].structure != NULL; strucs++){
    printf("%s %6.2f\n",sol[strucs].structure,sol[strucs].energy);
    gsl_histogram_increment(h,sol[strucs].energy);
    // TODO: so_structs[bin] += 1; ueber eigenes s histogram
    free(sol[strucs].structure);
  }
  free(sol);
}


