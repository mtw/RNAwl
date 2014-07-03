/*
  wl_options.c Command-line parsing for Wang-Landau sampling
  Last changed Time-stamp: <2014-07-03 18:09:32 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wl_options.h"
#include "wl_cmdline.h"
#include "ViennaRNA/utils.h"


static void ini_globals(void);
static void set_parameters(void);
static void display_settings(void);
static void to_basename(char *arg);
static void parse_infile(FILE *fp);

static struct gengetopt_args_info args_info;

/* ==== */
void 
process_commandline (int argc, char *argv[])
{
  ini_globals();
  
  if (cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  }
  set_parameters();
 
  if (args_info.inputs_num){
    char *infile =NULL;
    infile = strdup(args_info.inputs[0]);
    to_basename(args_info.inputs[0]);
    wanglandau_opt.INFILE = fopen(infile, "r");
    free(infile);
  }
  else
    wanglandau_opt.INFILE = stdin;
  
  parse_infile(wanglandau_opt.INFILE);
  
    if (args_info.inputs_num){
    fclose(wanglandau_opt.INFILE);
  }
  
}

/* ==== */
static void
ini_globals(void)
{
  wanglandau_opt.INFILE     = NULL;
  wanglandau_opt.bins       = 100;
  wanglandau_opt.ffinal     = 1e-7;
  wanglandau_opt.flat       = 0.8;
  wanglandau_opt.res        = 0.5;         /* kcal/mol */
  wanglandau_opt.steps      = 100000;
  wanglandau_opt.seed       = 123456789;
  wanglandau_opt.seed_given = 0;
  wanglandau_opt.T          = 37.0;
  wanglandau_opt.erange     = -1;
  wanglandau_opt.norm       = 1;
  wanglandau_opt.max        = 99999999999999.;
  wanglandau_opt.max_given  = 0;
  wanglandau_opt.verbose    = 0;
  wanglandau_opt.debug      = 0;
}

/* ==== */
/* process command line options */
static void
set_parameters(void)
{
  if (args_info.Temp_given){
    if( (wanglandau_opt.T = args_info.Temp_arg) < -273.15 ){
      fprintf(stderr, "Value of --Temp must be > -273.15\n");
      exit (EXIT_FAILURE);
    }
  }

   if(args_info.resolution_given){
     wanglandau_opt.res_given = 1;
    if( (wanglandau_opt.res = args_info.resolution_arg) < 0.1 ){
      fprintf(stderr, "Value of --resolution must be >= 0.1 \n");
      exit (EXIT_FAILURE);
    }
  }
   
  if (args_info.bins_given){
    if( (wanglandau_opt.bins = args_info.bins_arg) <= 0 ){
      fprintf(stderr, "Value of --bins must be > 0\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.steps_given){
    if( (wanglandau_opt.steps = args_info.steps_arg) <= 0 ){
      fprintf(stderr, "Value of --steps must be > 0\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.seed_given){
    wanglandau_opt.seed_given = 1;
    if( (wanglandau_opt.seed = args_info.seed_arg) <= 0 ){
      fprintf(stderr, "Value of --seed must be > 0\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if(args_info.mod_given){
    if( (wanglandau_opt.ffinal = args_info.mod_arg) < 1e-15 ){
      fprintf(stderr, "Value of -f must be > 1e-15 \n");
      exit (EXIT_FAILURE);
    }
  }
  
  if(args_info.flat_given){
    if( (wanglandau_opt.flat = args_info.flat_arg) <= 0.1 ){
      fprintf(stderr, "Value of --flat must be >= 0.1 \n");
      exit (EXIT_FAILURE);
    }
  }
  
  if(args_info.max_given){
    wanglandau_opt.max = args_info.max_arg;
    wanglandau_opt.max_given = 1;
  }

  if (args_info.norm_given){
    if( (wanglandau_opt.norm = args_info.norm_arg) < 2){
      fprintf(stderr, "Value of --norm must be >= 2 \n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.verbose_given){wanglandau_opt.verbose = 1;}
  if (args_info.debug_given){wanglandau_opt.debug = 1;}
  
  if (args_info.info_given){
    display_settings();
    exit(EXIT_SUCCESS);
  }
  
}

/* ==== */
static void
display_settings(void)
{
  fprintf(stderr, "Settings:\n");
  fprintf(stderr,
	  "--bins    = %d\n"
	  "--max     = %g\n"
	  "--ffinal  = %g\n"
	  "--flat    = %g\n"
	  "--norm    = %d\n"
	  "--res     = %g\n"
	  "--seed    = %lu\n"
	  "--steps   = %lu\n"
	  "--Temp    = %4.2f\n"
	  "--verbose = %i\n"
	  "--debug   = %i\n",
	  wanglandau_opt.bins,
	  wanglandau_opt.max,
	  wanglandau_opt.ffinal,
	  wanglandau_opt.flat,
	  wanglandau_opt.norm,
	  wanglandau_opt.res,
	  wanglandau_opt.seed,
	  wanglandau_opt.steps,
	  wanglandau_opt.T,
	  wanglandau_opt.verbose,
	  wanglandau_opt.debug);
}


/* ==== */
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
  wanglandau_opt.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  wanglandau_opt.structure = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(wanglandau_opt.sequence != NULL); assert(wanglandau_opt.structure != NULL);
  sscanf(line, "%s", wanglandau_opt.sequence);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", wanglandau_opt.structure);
  free (line);
  wanglandau_opt.len = strlen(wanglandau_opt.sequence);
}

/* ==== */
static void
to_basename(char *arg)
{
  int len;
  char *s=NULL, *t=NULL;
  
  s = strdup(arg);
  len = strlen(s);
  t = rindex(s, '/');
  if (t != NULL) memmove(s, t+1, (len-(t-s))*sizeof(char));
  t = NULL; t = index(s, '.');
  if (t != NULL) *t = '\0';
  wanglandau_opt.basename = strdup(s);
  free(s);
}


/* ==== */
void
dealloc_gengetopt(void)
{
  cmdline_parser_free (&args_info);
}

/* End of file */
