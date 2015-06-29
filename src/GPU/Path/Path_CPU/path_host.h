/*
 * path_host.h
 *
 *  Created on: Dec 6, 2014
 *      Author: yxc
 */

#ifndef PATH_HOST_H_
#define PATH_HOST_H_

#include "newton_host.h"
#include "predictor_host.h"

bool path_tracker(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, Parameter path_parameter,\
		          double& timeSec_Predict, double& timeSec_Eval, double& timeSec_MGS,\
		          int reverse = 0);

#endif /* PATH_HOST_H_ */
