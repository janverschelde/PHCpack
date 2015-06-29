/*
 * predict_test.h
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#ifndef PREDICT_TEST_H_
#define PREDICT_TEST_H_

#include "eval_host.h"
#include "path_gpu.h"
#include "err_check.h"

#include "predictor_host.h"

T1 predict_test(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, CT t);


#endif /* PREDICT_TEST_H_ */
