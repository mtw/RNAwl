/*  Last changed Time-stamp: <2014-05-10 00:27:38 mtw> */

#ifndef WL_RNA_H
#define WL_RNA_H

#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/subopt.h>


/* RNA-related */
void initialize_RNA(const char *);
void pre_process_RNA(void);
void post_process_RNA(void);

#endif
