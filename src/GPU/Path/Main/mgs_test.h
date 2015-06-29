/*
 * mgs_test.h
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#ifndef MGS_TEST_H_
#define MGS_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

#include "mgs_host.h"

T1 mgs_test(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom);

T1 mgs_test_large(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom);

T1 mgs_test_any(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, int device_option = 0, int n_path=-1);


#endif /* MGS_TEST_H_ */
