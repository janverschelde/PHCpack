/*
 * err_check.cpp
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#include "err_check.h"

T1 err_check_workspace(const CT* workspace1, const CT* workspace2, int n_workspace_size, int n_err_print) {
	T1 err(0);
	T1 err_bond(1E-10);
	int pr = 2 * sizeof(T1);
	std::cout.precision(pr);
	int err_count = 0;

	/*for(int i=0; i<n_workspace_size; i++) {
		std::cout << "d" << i << " = " << workspace1[i] \
		         << "      " << workspace2[i];
	}
	std::cout << "---------err------------------------------------" << std::endl;*/
	for(int i=0; i<n_workspace_size; i++) {
		bool err_print = 0;
		T1 tmp_err;
		if(workspace1[i].real == workspace1[i].real && workspace2[i].real == workspace2[i].real){
			tmp_err = abs(workspace1[i].real - workspace2[i].real);
		}
		else{
			tmp_err = max(abs(workspace1[i].real), abs(workspace2[i].real));
		}
		if(tmp_err > err_bond || workspace1[i].real != workspace1[i].real || workspace2[i].real != workspace2[i].real) {
			err_count++;
			if(err_count < n_err_print) {
				std::cout << "d" << i << " = " << workspace1[i]
				<< "     " << workspace2[i];
				err_print = 1;
			}
		}
		if(tmp_err > err) {
			err = tmp_err;
		}

		if(workspace1[i].imag == workspace1[i].imag && workspace2[i].imag == workspace2[i].imag){
			tmp_err = abs(workspace1[i].imag - workspace2[i].imag);
		}
		else{
			tmp_err = max(abs(workspace1[i].imag), abs(workspace2[i].imag));
		}
		if(tmp_err > err_bond || workspace1[i].imag != workspace1[i].imag || workspace2[i].imag != workspace2[i].imag) {
			err_count++;
			if((err_count < n_err_print) && (err_print == 0)) {
				std::cout << "d" << i << " = " << workspace1[i]
				<< "      " << workspace2[i];
			}
		}

		if(tmp_err > err) {
			err = tmp_err;
		}

	}
	if(err_count > 0) {
		std::cout << "n_err = " << err_count << std::endl;
	}
	//std::cout << "err = " << std::scientific << err << std::endl;
	return err;
}

T1 err_check_workspace_matrix(const CT* workspace1, const CT* workspace2, int n_rows, int n_cols) {
	T1 err(0);
	T1 err_bond(1E-10);
	int err_count = 0;
	int n_err_print = 20;
	int i=0;
	for(int col=0; col<n_cols; col++){
		for(int row=0; row<n_rows; row++) {
			bool err_print = 0;
			T1 tmp_err;
			if(workspace1[i].real == workspace1[i].real && workspace2[i].real == workspace2[i].real){
				tmp_err = abs(workspace1[i].real - workspace2[i].real);
			}
			else{
				tmp_err = max(abs(workspace1[i].real), abs(workspace2[i].real));
			}
			if(tmp_err > err_bond || workspace1[i].real != workspace1[i].real || workspace2[i].real != workspace2[i].real) {
				err_count++;
				if(err_count < n_err_print) {
					std::cout << "col="<< col << " "<< "row="<< row << " " << workspace1[i]\
					                               << "              " << workspace2[i];
					err_print = 1;
				}
			}
			if(tmp_err > err) {
				err = tmp_err;
			}

			if(workspace1[i].imag == workspace1[i].imag && workspace2[i].imag == workspace2[i].imag){
				tmp_err = abs(workspace1[i].imag - workspace2[i].imag);
			}
			else{
				tmp_err = max(abs(workspace1[i].imag), abs(workspace2[i].imag));
			}
			if(tmp_err > err_bond || workspace1[i].imag != workspace1[i].imag || workspace2[i].imag != workspace2[i].imag) {
				err_count++;
				if((err_count < n_err_print) && (err_print == 0)) {
					std::cout << "col="<< col << " "<< "row="<< row << " " << workspace1[i]\
					                               << "              " << workspace2[i];
				}
			}

			if(tmp_err > err) {
				err = tmp_err;
			}
			i++;
		}
	}
	if(err_count > 0) {
		std::cout << "n_err = " << err_count << std::endl;
	}
	return err;
}

void err_check_class_workspace(CT** deri_val, CT* f_val, CT* matrix, int n_eq, int dim) {
	T1 err(0);
	T1 err_bond(1E-10);
	for(int i=0; i<n_eq; i++) {
		for(int j=0; j<dim; j++) {
			T1 tmp_err = abs(deri_val[i][j].real - matrix[j*n_eq + i].real);
			//std::cout << i << " " << j  << " " << tmp_err << std::endl;
			if(tmp_err > err_bond) {
				std::cout << "d" << i << " " << j << " = " << deri_val[i][j]
				<< "       " << matrix[j*n_eq + i];
			}
			if(tmp_err > err) {
				err = tmp_err;
			}
			tmp_err = abs(deri_val[i][j].imag - matrix[j*n_eq + i].imag);
			if(tmp_err > err_bond) {
				err = tmp_err;
				std::cout << "d" << i << " " << j << " = " << deri_val[i][j]
				<< "       " << matrix[j*n_eq + i];
			}
			if(tmp_err > err) {
				err = tmp_err;
			}

		}
		T1 tmp_err = abs(f_val[i].real - matrix[dim*n_eq + i].real);
		if(tmp_err > err_bond) {
			std::cout << "f" << i << " = " << f_val[i]
			<< "     " << matrix[dim*n_eq + i];
		}
		if(tmp_err > err) {
			err = tmp_err;
		}
		tmp_err = f_val[i].imag - matrix[dim*n_eq + i].imag;
		if(tmp_err > err_bond) {
			err = tmp_err;
			std::cout << "f" << i << " = " << f_val[i]
			<< "    " << matrix[dim*n_eq + i];
		}
		if(tmp_err > err) {
			err = tmp_err;
		}
	}
	//std::cout << "err = " << std::scientific << err << std::endl;
}

T1 err_check_r(CT** CPU_R, CT* GPU_R, int dim, int right_hand_side) {
	T1 err(0);
	T1 err_bond(1E-10);
	int err_count = 0;
	int n_err_print = 10;

	int rows;
	int cols;

	if(right_hand_side == 1){
		rows = dim+1;
		cols = dim+1;
	}
	else{
		rows = dim;
		cols = dim;
	}

	for(int j=0; j<cols; j++){
		for(int i=0; i<=j; i++){
			if(i==dim && j == dim && right_hand_side == 1){
				break;
			}
			int tmp_idx =(dim+1)*(dim+2)/2 -(j+2)*(j+1)/2 + i;
			CT tmp_cpu = CPU_R[i][j];
			CT tmp_gpu = GPU_R[tmp_idx];

			//std::cout << i << " " << j << " " << GPU_R[tmp_idx];
			//std::cout  << i << " " << j << " " << CPU_R[i][j];

			T1 tmp_err;
			if(tmp_cpu.real == tmp_cpu.real && tmp_gpu.real == tmp_gpu.real){
				tmp_err = abs(tmp_cpu.real - tmp_gpu.real);
			}
			else{
				tmp_err = max(abs(tmp_cpu.real), abs(tmp_gpu.real));
			}
			if(tmp_err > err_bond) {
				err_count++;
				if(err_count < n_err_print) {
					std::cout << "d" << i << " " << j << " = " << tmp_cpu
					<< "     " << tmp_gpu;
				}
			}
			if(tmp_err > err) {
				err = tmp_err;
			}

			if(tmp_cpu.imag == tmp_cpu.imag && tmp_gpu.imag == tmp_gpu.imag){
				tmp_err = abs(tmp_cpu.imag - tmp_gpu.imag);
			}
			else{
				tmp_err = max(abs(tmp_cpu.imag), abs(tmp_gpu.imag));
			}
			if(tmp_err > err_bond) {
				err_count++;
				if(err_count < n_err_print) {
					std::cout << "d" << i << " " << j << " = " << tmp_cpu
					<< "     " << tmp_gpu;
				}
			}

			if(tmp_err > err) {
				err = tmp_err;
			}
		}
	}

	/*for(int i=0; i<n_workspace_size; i++) {
		T1 tmp_err;
		if(workspace1[i].real == workspace1[i].real && workspace2[i].real == workspace2[i].real){
			tmp_err = abs(workspace1[i].real - workspace2[i].real);
		}
		else{
			tmp_err = max(abs(workspace1[i].real), abs(workspace2[i].real));
		}
		if(tmp_err > err_bond) {
			err_count++;
			if(err_count < n_err_print) {
				std::cout << "d" << i << " = " << workspace1[i]
				<< "     " << workspace2[i];
			}
		}
		if(tmp_err > err) {
			err = tmp_err;
		}

		if(workspace1[i].imag == workspace1[i].imag && workspace2[i].imag == workspace2[i].imag){
			tmp_err = abs(workspace1[i].imag - workspace2[i].imag);
		}
		else{
			tmp_err = max(abs(workspace1[i].imag), abs(workspace2[i].imag));
		}
		if(tmp_err > err_bond) {
			err_count++;
			if(err_count < n_err_print) {
				std::cout << "d" << i << " = " << workspace1[i]
				<< "     " << workspace2[i];
			}
		}

		if(tmp_err > err) {
			err = tmp_err;
		}
	}*/

	if(err_count > 0) {
		std::cout << "n_err = " << err_count << std::endl;
	}
	return err;
}
