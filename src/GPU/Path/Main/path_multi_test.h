/*
 * multi_path_test.h
 *
 *  Created on: Feb 12, 2015
 *      Author: yxc
 */

#ifndef PATH_MULTI_TEST_H_
#define PATH_MULTI_TEST_H_

#include "path_test.h"

bool path_multiple_test(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, Parameter path_parameter,\
		PolySolSet& sol_set, PolySys& Target_Sys, int device_option, int n_path=-1);


#endif /* PATH_MULTI_TEST_H_ */
