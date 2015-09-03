/* path_gpu.h, created on Dec 3, 2014 by yxc with edits by jv */

#ifndef EVAL_KERNEL_H_
#define EVAL_KERNEL_H_

#include <iostream>

#include "DefineType_Host.h"
#include "eval_host.h"
#include "parameter.h"

int GPU_Eval
 ( const CPUInstHom& hom, CT* cpu_sol0, CT* cpu_t, CT**& gpu_workspace_all,
   CT**& gpu_matrix_all, int n_path = 1, int* x_t_idx=NULL,
   int n_predictor=1 );

int GPU_Predict
 ( const CPUInstHom& hom, CT*& x_gpu, int n_predictor, CT cpu_t,
   int n_path = 1 );

int GPU_MGS_Mult
 ( const CPUInstHom& hom, CT**& sol_gpu, CT**& matrix_gpu_q,
   CT**& matrix_gpu_r, CT* V, int n_path=1 );

int GPU_MGS
 ( const CPUInstHom& hom, CT*& sol_gpu, CT*& matrix_gpu_q, CT*& matrix_gpu_r,
   CT* V );

bool GPU_Newton
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT cpu_t,
   CT*& x_new, int n_path=1, int verbose=0 );

bool GPU_Path
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT cpu_t,
   CT*& x_gpu, int n_path=1, int inverse=0, int verbose=0 );

bool* GPU_Path_mult
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT* cpu_t,
   CT**& x_gpu, int n_path=1, int inverse=0, int verbose=0 );

#endif /* EVAL_KERNEL_H_ */
