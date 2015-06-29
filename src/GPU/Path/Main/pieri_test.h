/*
 * pieri_test.h
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#ifndef PIERI_TEST_H_
#define PIERI_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

//#include "path_host.h"
#include "init_test.h"
#include "path_test.h"

void Pieri_Test(int dim_start, int dim_end, Parameter path_parameter, int sys_type, int device_option);

#endif /* PIERI_TEST_H_ */
