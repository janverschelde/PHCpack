/*
 * err_check.h
 *
 *  Created on: Feb 8, 2015
 *      Author: yxc
 */

#ifndef ERR_CHECK_H_
#define ERR_CHECK_H_

#include <iostream>
#include "DefineType.h"

T1 err_check_r(CT** CPU_R, CT* GPU_R, int dim, int right_hand_side=1);

T1 err_check_workspace(const CT* workspace1, const CT* workspace2, int n_workspace_size, int n_err_print = 20);

T1 err_check_workspace_matrix(const CT* workspace1, const CT* workspace2, int n_rows, int n_cols);

void err_check_class_workspace(CT** deri_val, CT* f_val, CT* matrix, int n_eq, int dim);

#endif /* ERR_CHECK_H_ */
