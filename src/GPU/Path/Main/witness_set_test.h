/*
 * witness_set_test.h
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#ifndef WITNESS_SET_TEST_H_
#define WITNESS_SET_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"

#include "path_test.h"

#include <vector>

int witness_set_test(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, Parameter path_parameter, CT* sol0, CT t, int device_option = 1);



#endif /* WITNESS_SET_TEST_H_ */
