/*
  wl_options.c Command-line parsing for Wang-Landau sampling
  Last changed Time-stamp: <2014-06-26 17:23:38 mtw>
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
static void warn (FILE *hdl, char *fmt, ...);

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
}

/* ==== */
static void
ini_globals(void)
{
  wanglandau_opt.INFILE  = NULL;
  wanglandau_opt.bins    = 30;
  wanglandau_opt.ffinal  = 1e-7;
  wanglandau_opt.flat    = 0.8;
  wanglandau_opt.steps   = 1;
  wanglandau_opt.T       = 37.0;
  wanglandau_opt.erange  = -1;
  wanglandau_opt.norm    = 1;
  wanglandau_opt.verbose = 0;
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
  
  if(args_info.flat_given){
    if( (wanglandau_opt.flat = args_info.flat_arg) <= 0.1 ){
      fprintf(stderr, "Value of --flat must be >= 0.1 \n");
      exit (EXIT_FAILURE);
    }
  }

  
  
}

/* ==== */
static void
display_settings(void)
{
  fprintf(stderr, "Settings:\n");
  fprintf(stderr,
	  "--bins    = %d\n"
	  "--Temp    = %4.2f\n",
	  wanglandau_opt.bins,
	  wanglandau_opt.T
	  );
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
static void
warn (FILE *hdl, char *fmt, ...)
{
  va_list args;
  if ( hdl == NULL ) hdl = stderr;
  va_start(args, fmt);
  fprintf(hdl, "WARNING: ");
  vfprintf(hdl, fmt, args);
  fprintf(hdl,"\n");
  fflush(hdl);
  va_end(args);
}

/* ==== */
void
dealloc_gengetopt(void)
{
  cmdline_parser_free (&args_info);
}

/* End of file */
