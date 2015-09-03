/* Parameter_Header.h created on 1 Feb 2015 by yxc with edits by jv */

#ifndef PARAMETER_HEADER_H_
#define PARAMETER_HEADER_H_

#define ERR 1E-10

#define maxrounds 128

#define max_array_size 2000

#define shmemsize 512

#define MON_EVAL_METHOD 2
// 0 : reverse mode
// 1 : reverse mode with aligned memory for instructions
// 2 : tree mode
// 3 : for multiple evaluation (chosen when #paths > 1)

#define BS_QR 256

#define BS_QR_Back 256

#define BS_Mon_Align 64

// QR Reduce Parameters
#define matrix_block_row 32
#define matrix_block_pivot_col 4
#define matrix_block_reduce_col 4

// Parameters
#define N_PREDICTOR           4

#define MAX_STEP              400
#define MAX_DELTA_T           1E-1
#define MAX_DELTA_T_END       1E-2
#define MIN_DELTA_T           1E-7

#define MAX_IT                3
#define ERR_MIN_ROUND_OFF     1E-9

#define MAX_IT_REFINE         		5
#define ERR_MIN_ROUND_OFF_REFINE    1E-11

#define ERR_MAX_RES           1E-2
#define ERR_MAX_DELTA_X       1E-1
#define ERR_MAX_FIRST_DELTA_X 1E-2

#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7

#endif /* PARAMETER_HEADER_H_ */
