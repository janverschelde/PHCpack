/* Parameter_Header.h created on 1 Feb 2015 by yxc with edits by jv */

#ifndef PARAMETER_HEADER_QD_H_
#define PARAMETER_HEADER_QD_H_

#define ERR 1E-55

#define maxrounds 128

#define max_array_size 2000

#define shmemsize 128

#define MON_EVAL_METHOD 1
// 0 : reverse mode
// 1 : reverse mode with aligned memory for instructions
// 2 : tree mode
// 3 : for multiple evaluation (chosen when #paths > 1)

#define BS_QR 64

#define BS_QR_Back 64

#define BS_Mon_Align 128

#define matrix_block_row 32
#define matrix_block_pivot_col 1
#define matrix_block_reduce_col 1

// Parameters
#define N_PREDICTOR 4

#define MAX_STEP              2000
#define MAX_DELTA_T           1E-1
#define MIN_DELTA_T           1E-20
#define MAX_DELTA_T_END       1E-2

#define MAX_IT                7
#define ERR_MAX_RES           1E-1
#define ERR_MAX_DELTA_X       1E-1
#define ERR_MAX_FIRST_DELTA_X 1E-2
#define ERR_MIN_ROUND_OFF     1E-30

#define MAX_IT_REFINE         		9
#define ERR_MIN_ROUND_OFF_REFINE    1E-40

#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7

#endif /* PARAMETER_HEADER_QD_H_ */
