/*
  Last changed Time-stamp: <2014-05-05 11:24:19 mtw>
  Literature:
  Landau, PD and Tsai, S-H and Exler, M (2004) Am. J. Phys. 72:(10) 1294-1302
  A new approach to Monte Carlo simulations in statistical physics:
  Wang-Landau sampling
*/

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include "globals.h"

int
main (int argc, char **argv)
{
  if ( signal (SIGUSR1, &sighandler) == SIG_ERR)
    fprintf(stderr,"Couldn't register signal handler\n");

  process_commandline(argc,argv);
  wanglandau();
  wanglandau_free_memory();

  return (EXIT_SUCCESS);
}

/* End of file */
