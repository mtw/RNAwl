#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wl_options.h"
#include "ViennaRNA/utils.h"


static void ini_globals(void);
//static void set_parameters(void);
//static void display_settings(void);
//static void to_basename(char *arg);
static void parse_infile(FILE *fp);
//static void  warn (FILE *hdl, char *fmt, ...);


/* ==== */
void 
process_commandline (int argc, char *argv[])
{
  ini_globals();
  /*
  if (cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  }
  set_parameters();
 
  if (args_info.inputs_num){
    char *infile =NULL;
    infile = strdup(args_info.inputs[0]);
    to_basename(args_info.inputs[0]);
    wl_opt.INFILE = fopen(infile, "r");
    free(infile);
  }
  else
    wl_opt.INFILE = stdin;
  */
  parse_infile(wanglandau_opt.INFILE);
}

/* ==== */
static void
ini_globals(void)
{
  wanglandau_opt.INFILE = stdin;
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
