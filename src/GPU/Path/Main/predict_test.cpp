/*
 * predict_test.cpp
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#include "predict_test.h"

T1 predict_test(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, CT t) {
	std::cout << "--------- Predict Test ----------" << std::endl;
	workspace_cpu.init_x_t_predict_test();
	workspace_cpu.update_t_value(t);
	predictor_newton(workspace_cpu.x_array, workspace_cpu.t_array,
	workspace_cpu.x_t_idx, cpu_inst_hom.n_predictor, cpu_inst_hom.dim);
	CT* x_gpu;
	GPU_Predict(cpu_inst_hom, x_gpu, cpu_inst_hom.n_predictor, t);

	std::cout << "--------- GPU Predictor Error----------" << std::endl;
	T1 err = err_check_workspace(workspace_cpu.x, x_gpu, cpu_inst_hom.dim);
	std::cout << "x_cpu[0] = " << workspace_cpu.x[0];
	std::cout << "x_gpu[0] = " << x_gpu[0];
	free(x_gpu);
	return err;
}


