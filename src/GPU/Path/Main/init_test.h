/* init_test.h, created on Feb 8, 2015 by yxc with edits by jv */

#ifndef INIT_TEST_H_
#define INIT_TEST_H_

#include "families.h"
#include "eval_host.h"
#include "polysol.h"

bool init_homotopy_test
 ( PolySys& Target_Sys, PolySys& Start_Sys, PolySolSet& sol_set,
   int dim, int& n_eq, CT*& sol0,
   CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu, int test, 
   int n_predictor, string start_sys_filename, string target_sys_filename );

void init_cpu_inst_workspace
 ( PolySys& Target_Sys, PolySys& Start_Sys, int dim, int n_eq,
   int n_predictor, CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu,
   int test );

bool read_homotopy_file
 ( PolySys& Target_Sys, PolySys& Start_Sys, int dim, int& n_eq,
   string Start_Sys_filename, string Target_Sys_file_name,
   PolySolSet* sol_set=NULL );

bool init_homotopy_test_ada
 ( PolySys& Target_Sys, PolySys& Start_Sys, 
   CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu, int n_predictor );

#endif /* INIT_TEST_H_ */
