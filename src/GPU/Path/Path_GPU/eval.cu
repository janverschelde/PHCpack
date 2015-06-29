#if(path_precision == 0)
#include "eval_mon_d.cu"
#elif(path_precision == 1)
#include "eval_mon_dd.cu"
#else
#include "eval_mon_qd.cu"
#endif

#include "eval_sum.cu"
#include "eval_coef.cu"
#include "eval_mult.cu"
#include "path_init.cu"

void eval(GPUWorkspace& workspace, const GPUInst& inst) {

	//std::cout << "workspace.n_path_continuous = " << workspace.n_path_continuous << std::endl;
	/*dim3 coef_grid = get_grid(inst.n_coef, inst.coef_BS, workspace.n_path_continuous);
	eval_coef_kernel<<<coef_grid, inst.coef_BS>>>(workspace.coef,\
			inst.coef, inst.n_coef, workspace.t_array, workspace.one_minor_t,\
			workspace.workspace_size, workspace.x_t_idx_mult, workspace.path_idx);*/

	//std::cout << "inst.n_coef = " << inst.n_coef << std::endl;
	eval_coef_kernel<<<inst.coef_grid, inst.coef_BS>>>(workspace.coef,\
			inst.coef, inst.n_coef, workspace.t, workspace.one_minor_t,\
			workspace.workspace_size);

	/*CT** gpu_workspace_all = new CT*[workspace.n_path];
	CT** gpu_matrix_all = new CT*[workspace.n_path];
	for(int path_idx=0; path_idx<workspace.n_path; path_idx++){
		gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
		gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
		std::cout << "path_idx = " << path_idx << std::endl;
		for(int coef_idx=0; coef_idx<10; coef_idx++){
			std::cout << coef_idx << " " <<  gpu_workspace_all[path_idx][coef_idx];
		}
	}*/

	eval_mon(workspace, inst);

	eval_sum(workspace, inst);
}

int GPU_Eval(const CPUInstHom& hom, CT* cpu_sol0, CT* cpu_t, \
		CT**& gpu_workspace_all, CT**& gpu_matrix_all, int n_path, int* x_t_idx, int n_predictor) {

	cout << "GPU Eval" << endl;
	// CUDA configuration
	cuda_set();

	GPUInst inst(hom, n_path);
	GPUWorkspace workspace(inst.mon_pos_size, inst.n_coef, inst.n_constant, inst.n_eq, inst.dim, n_predictor, inst.alpha, n_path);

	dim3 update_t_grid = get_grid(n_path,inst.predict_BS,1);
	std::cout << "n_path = " << n_path << std::endl;
	std::cout << "inst.n_sum_zero = " << inst.n_sum_zero << std::endl;
	path_mult_init_kernel<<<update_t_grid, inst.predict_BS>>>(workspace.t_mult, workspace.t_last_mult, \
			workspace.delta_t_mult, workspace.t_array, workspace.one_minor_t, workspace.alpha_gpu, \
			workspace.path_success, workspace.newton_success, workspace.end_range, workspace.n_success, workspace.max_f_val_last_gpu, \
			workspace.n_point_mult, workspace.x_t_idx_mult, workspace.workspace_size, workspace.path_idx, n_path, \
			workspace.matrix, 0, inst.sum_zeros, workspace.n_predictor);

	if(x_t_idx!=NULL){
		workspace.update_x_t_idx_all(x_t_idx);
		workspace.update_x_t_value_array(cpu_sol0, cpu_t, x_t_idx);
	}
	else{
		workspace.update_x_t_value(cpu_sol0, *cpu_t);
	}
	workspace.update_x_value_mult2(cpu_sol0, x_t_idx);
	workspace.update_t_value_mult2(cpu_t, x_t_idx);
	//workspace.init_workspace_eq(inst.n_pos_total_eq, workspace.n_path);

	bool mult = false;
	if(MON_EVAL_METHOD == 3){
		mult = true;
	}
	if(mult){
		eval_mult(workspace, inst);
		//std::cout << "n_workspace = " << inst.n_workspace << std::endl;
		gpu_workspace_all = new CT*[n_path];
		gpu_matrix_all = new CT*[n_path];
		CT* gpu_matrix_mult = workspace.get_matrix_mult();
		CT* gpu_workspace_mult = workspace.get_workspace();
		for(int path_idx=0; path_idx<n_path; path_idx++){
			gpu_matrix_all[path_idx] = new CT[workspace.n_matrix];
			for(int i=0; i<workspace.n_matrix; i++){
				gpu_matrix_all[path_idx][i] = gpu_matrix_mult[n_path*i+path_idx];
			}
			int n = workspace.n_coef+workspace.mon_pos_size;
			gpu_workspace_all[path_idx] = new CT[n];
			for(int i=0; i<n; i++){
				gpu_workspace_all[path_idx][i] = gpu_workspace_mult[n_path*i+path_idx];
			}
		}
	}
	else{
		eval(workspace, inst);
		//std::cout << "n_workspace = " << inst.n_workspace << std::endl;
		gpu_workspace_all = new CT*[n_path];
		gpu_matrix_all = new CT*[n_path];
		for(int path_idx=0; path_idx<n_path; path_idx++){
			gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
			gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
		}
	}

	cudaDeviceReset();
	return 0;
}
