/* multi_path_test.h, created on Feb 12, 2015 by yxc with edits by jv */

#ifndef PATH_MULTI_TEST_H_
#define PATH_MULTI_TEST_H_

#include "path_test.h"
#include "polysol.h"

bool path_multiple_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, PolySolSet& sol_set, PolySys& Target_Sys,
   int device_option, int n_path=-1 );

#endif /* PATH_MULTI_TEST_H_ */
