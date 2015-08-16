/* path_test.h, created on Feb 8, 2015 by yxc with edits by jv */

#ifndef PATH_TEST_H_
#define PATH_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

#include "path_host.h"

bool path_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT*& sol_new, CT t,
   PolySys& Target_Sys, int device_option = 0 );

T1 path_test_reverse
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT t, int device_option = 1 );

#endif /* PATH_TEST_H_ */
