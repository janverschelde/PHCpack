/* newton_test.h, created on Feb 8, 2015 by yxc with edits by jv */

#ifndef NEWTON_TEST_H_
#define NEWTON_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

#include "newton_host.h"

T1 newton_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT t );
/*
 * DESCRIPTION :
 *   Tests the calling of Newton's method. */

#endif /* NEWTON_TEST_H_ */
