/*
 * newton_host.cpp
 *
 *  Created on: Dec 6, 2014
 *      Author: yxc
 */

#include "newton_host.h"

double norm1(double real, double imag){
	return abs(real) + abs(imag);
}

double CPU_normalize(complexH<double>* sol, int dim){
	double max_delta = norm1(sol[0].real,sol[0].imag);
	for(int k=1; k<dim; k++){
		double tmp_delta = norm1(sol[k].real,sol[k].imag);
		if(tmp_delta>max_delta){
			max_delta = tmp_delta;
		}
	}
	return max_delta;
}

double CPU_normalize(complexH<dd_real>* sol, int dim){
	double max_delta = norm1(sol[0].real.x[0],sol[0].imag.x[0]);
	for(int k=1; k<dim; k++){
		double tmp_delta = norm1(sol[k].real.x[0],sol[k].imag.x[0]);
		if(tmp_delta>max_delta){
			max_delta = tmp_delta;
		}
	}
	return max_delta;
}

double CPU_normalize(complexH<qd_real>* sol, int dim){
	double max_delta = norm1(sol[0].real.x[0],sol[0].imag.x[0]);
	for(int k=1; k<dim; k++){
		double tmp_delta = norm1(sol[k].real.x[0],sol[k].imag.x[0]);
		if(tmp_delta>max_delta){
			max_delta = tmp_delta;
		}
	}
	return max_delta;
}

bool CPU_Newton(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, Parameter path_parameter,
                double& timeSec_Eval, double& timeSec_MGS, int reverse){
    //cout << "Newton max_it = " << path_parameter.max_it << " err_max_delta_x = " << path_parameter.err_max_delta_x << endl;

    //clock_t begin = clock();


    bool Debug = false;
    //Debug = true;
    if(workspace_cpu.path_idx == 0){
    	//Debug = true;
    }

	bool Record = false;
	//Record = true;

    CT* x = workspace_cpu.x;
    CT t = *(workspace_cpu.t);
    CT** V = (workspace_cpu.V);
    CT** R = (workspace_cpu.R);
    CT* sol = (workspace_cpu.sol);
    int dim = cpu_inst_hom.dim;
    int n_eq = cpu_inst_hom.n_eq;

	// Parameters
    double err_round_off;
	if(t.real > 0.9){
		err_round_off = path_parameter.err_min_round_off_refine;
	}
	else{
		err_round_off = path_parameter.err_min_round_off;
	}

	//double last_max_delta_x = path_parameter.err_max_first_delta_x;
	//double last_max_f_val   = path_parameter.err_max_res;
	//double last_max_f_val   = 1E10;

    double max_x = CPU_normalize(x,dim);
    /*std::cout << "X value" << std::endl;
    for(int var_idx=0; var_idx<dim; var_idx++){
    	std::cout << var_idx << " " << x[var_idx];
    }
    std::cout << "max_x = " << max_x << std::endl;*/

    bool success = 0;

    //clock_t begin_eval = clock();
    //t = CT(0.0, 0.0);
	cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
	cpu_inst_hom.n_eval_CPU++;
    //clock_t end_eval = clock();
    //timeSec_Eval += (end_eval - begin_eval) / static_cast<double>( CLOCKS_PER_SEC );

    double max_f_val = CPU_normalize(V[dim],n_eq);
	double r_max_f_val = max_f_val/max_x;
	if(Debug){
		std::cout.precision(2);
		std::cout << scientific;
		std::cout << "          max_x : " << max_x << std::endl;
		std::cout << "   residual(a&r): " << max_f_val << " " << r_max_f_val << std::endl;
	}

    /*if(max_f_val != max_f_val){
    	return success;
    }*/

    /*if(r_max_f_val > 1E-2){
    	return success;
    }*/

    if(max_f_val < err_round_off || r_max_f_val < err_round_off){
    	success = 1;
    	return success;
    }

    double last_max_f_val = max_f_val;

    double max_delta_x;

    for(int i=0; i<path_parameter.max_it; i++){
	    //clock_t begin_mgs = clock();
    	/*if(Debug){
			std::cout << "matrix" << std::endl;
			int pr = 2 * sizeof(T1);
			std::cout.precision(pr);
			for(int var_idx=1; var_idx<=1; var_idx++){
				CT tmp(0.0,0.0);
				for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
					std::cout << var_idx << " " << eq_idx << " " \
							<< V[var_idx][eq_idx];
					tmp += V[var_idx][eq_idx];
				}
				std::cout << tmp;
			}
			std::cout << std::endl;
	    }*/
        CPU_mgs2qrls(V,R,sol,n_eq,dim+1);
    	cpu_inst_hom.n_mgs_CPU++;
	    //clock_t end_mgs = clock();
	    //timeSec_MGS += (end_mgs - begin_mgs) / static_cast<double>( CLOCKS_PER_SEC );

	    max_delta_x = CPU_normalize(sol,dim);
	    /*if(Debug){
	    	std::cout << "sol" << std::endl;
			for(int var_idx=0; var_idx<dim; var_idx++){
				std::cout << var_idx << " " << sol[var_idx];
			}
	    }*/

		double r_max_delta_x = max_delta_x/max_x;

    	if(Debug){
    		std::cout << " correction(a&r): " << max_delta_x  << " " << r_max_delta_x;
    	}

    	if(Record){
    		cpu_inst_hom.path_data.add_iteration(max_delta_x, r_max_delta_x);
    	}

        if(max_delta_x < err_round_off || r_max_delta_x < err_round_off){
        	success = 1;
        	break;
        }

        for(int k=0; k<dim; k++){
        	x[k] = x[k] - sol[k];
        }

 		/*std::cout << std::endl << "Correct X:" << std::endl;
	    for(int i=0; i<dim; i++){
	    	std::cout << i << " " << x[i];
	    }*/

	    max_x = CPU_normalize(x,dim);

	    //clock_t begin_eval = clock();
		cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
		cpu_inst_hom.n_eval_CPU++;
		//clock_t end_eval = clock();
		//timeSec_Eval += (end_eval - begin_eval) / static_cast<double>( CLOCKS_PER_SEC );

		max_f_val = CPU_normalize(V[dim],n_eq);
		r_max_f_val = max_f_val/max_x;
    	if(Debug){
    		std::cout << std::endl << "          max_x : " << max_x << std::endl;
    		std::cout << " residual(a&r): " << max_f_val << " " << r_max_f_val << std::endl;
    	}
    	if(Record){
    		cpu_inst_hom.path_data.update_iteration_res(max_f_val, r_max_f_val);
    	}

	    if(max_f_val > last_max_f_val){
	    	break;
	    }

		if(max_f_val < err_round_off || r_max_f_val < err_round_off){
			success = 1;
			break;
		}
        last_max_f_val = max_f_val;
    }

    /*if(success){
    	if(max_delta_x > path_parameter.err_max_delta_x){
    		//std::cout << "----------  Fail  -----------" << std::endl;
        	//std::cout << "Fail tolerance: " << last_max_delta_x << std::endl;
        	success = 0;
        }
    }*/

    //clock_t end = clock();
    //double timeSec = (end - begin) / static_cast<double>( CLOCKS_PER_SEC );
    //cout << "done: "<< timeSec << endl;

	if(Debug){
		std::cout << std::endl;
		std::cout << "Newton success = " << success << std::endl;
	}

	if(Record){
    	cpu_inst_hom.path_data.update_step_correct_pt(x);
	}
    return success;
}

bool CPU_Newton_Refine(Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, Parameter path_parameter,
                double& timeSec_Eval, double& timeSec_MGS, int reverse){
    //cout << "Newton max_it = " << path_parameter.max_it << " err_max_delta_x = " << path_parameter.err_max_delta_x << endl;

	// Parameters
    // eqs square to compare with normal square

    //clock_t begin = clock();

    CT* x = workspace_cpu.x;
    CT t = *(workspace_cpu.t);
    CT** V = (workspace_cpu.V);
    CT** R = (workspace_cpu.R);
    CT* sol = (workspace_cpu.sol);
    int dim = cpu_inst_hom.dim;
    int n_eq = cpu_inst_hom.n_eq; // to be removed

	double last_max_delta_x = path_parameter.err_max_first_delta_x;
	double last_max_f_val   = path_parameter.err_max_res;
    bool success = 1;

    bool Debug = false;
    //Debug = true;

    for(int i=0; i<path_parameter.max_it_refine; i++){
    	if(Debug){
    		std::cout << "  Iteration " << i << std::endl;
    	}
	    //clock_t begin_eval = clock();
    	cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
    	cpu_inst_hom.n_eval_CPU++;
	    //clock_t end_eval = clock();
	    //timeSec_Eval += (end_eval - begin_eval) / static_cast<double>( CLOCKS_PER_SEC );

	    double max_f_val = CPU_normalize(V[dim],n_eq);
    	if(Debug){
    		std::cout << "       max_f_value = " << max_f_val << std::endl;
    	}

        if(max_f_val > last_max_f_val || max_f_val != max_f_val){
        	success = 0;
        	break;
        }

        if(max_f_val < path_parameter.err_min_round_off_refine){
    		last_max_delta_x = 0;
        	break;
        }

	    //clock_t begin_mgs = clock();
        CPU_mgs2qrls(V,R,sol,n_eq,dim+1);
    	cpu_inst_hom.n_mgs_CPU++;
	    //clock_t end_mgs = clock();
	    //timeSec_MGS += (end_mgs - begin_mgs) / static_cast<double>( CLOCKS_PER_SEC );

	    double max_delta_x = CPU_normalize(sol,dim);

    	if(Debug){
    		std::cout << "       max_delta_x = " << max_delta_x << std::endl;
    	}

        if(max_delta_x > last_max_delta_x || max_delta_x != max_delta_x){
        	success = 0;
        	break;
        }

        if(max_delta_x < path_parameter.err_min_round_off){
        	last_max_delta_x = max_delta_x;
        	break;
        }

        last_max_delta_x   = max_delta_x;
        last_max_f_val = max_f_val;

        for(int k=0; k<dim; k++){
        	x[k] = x[k] - sol[k];
        }
    }

    if(success){
    	if(last_max_delta_x > path_parameter.err_max_delta_x){
    		//std::cout << "----------  Fail  -----------" << std::endl;
        	//std::cout << "Fail tolerance: " << last_max_delta_x << std::endl;
        	success = 0;
        }
    }

    cpu_inst_hom.max_delta_x = sqrt(last_max_delta_x);

    //clock_t end = clock();
    //double timeSec = (end - begin) / static_cast<double>( CLOCKS_PER_SEC );
    //cout << "done: "<< timeSec << endl;
    return success;
}


