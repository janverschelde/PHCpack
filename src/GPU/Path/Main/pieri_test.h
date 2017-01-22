/* The file "pieri_test.h" provide the main prototype to test path tracking
 * in a large dimensional space.  The paths are defined by a Pieri homotopy. 
 * The file was made by Xiangcheng Yu on 8 Febuary 2015. */

#ifndef PIERI_TEST_H_
#define PIERI_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

#include "init_test.h"
#include "path_test.h"

void Pieri_Test
 ( int dim_start, int dim_end, Parameter path_parameter, int sys_type,
   int device_option );

// Runs a sequence of path tracking jobs defined by Pieri homotopies.

#endif /* PIERI_TEST_H_ */
