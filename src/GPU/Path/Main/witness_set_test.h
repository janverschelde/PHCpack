/* The file "witness_set_test.h" was written by Xiangcheng Yu
 * on 8 February 2015. */

#ifndef WITNESS_SET_TEST_H_
#define WITNESS_SET_TEST_H_

#include <vector>
#include "eval_host.h"
#include "path_gpu.h"
#include "path_test.h"
#include "polysol.h"

int witness_set_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT t, int device_option = 1 );

#endif /* WITNESS_SET_TEST_H_ */
