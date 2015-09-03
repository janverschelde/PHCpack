/* Parameter_Header.h created on 1 Feb 2015 by yxc with edits by jv */

#ifndef PARAMETER_HEADER_DD_H_
#define PARAMETER_HEADER_DD_H_

#define ERR 1E-25

#define maxrounds 128

#define max_array_size 2000

#define shmemsize 256

#define MON_EVAL_METHOD 1
// 0 : reverse mode
// 1 : reverse mode with aligned memory for instructions
// 2 : tree mode
// 3 : for multiple evaluation (chosen when #paths > 1)

#define BS_QR 128

#define BS_QR_Back 128

#define BS_Mon_Align 128

#define matrix_block_row 64
#define matrix_block_pivot_col 1
#define matrix_block_reduce_col 2

// Parameters
#define N_PREDICTOR           4

#define MAX_STEP              1000
#define MAX_DELTA_T           1E-1
#define MAX_DELTA_T_END       1E-2
#define MIN_DELTA_T           1E-20

#define MAX_IT                5
#define ERR_MIN_ROUND_OFF     1E-15

#define MAX_IT_REFINE         		7
#define ERR_MIN_ROUND_OFF_REFINE    1E-25 //game8 -30 doesnt work

#define ERR_MAX_RES           1E-1
#define ERR_MAX_DELTA_X       1E-1
#define ERR_MAX_FIRST_DELTA_X 1E-2

#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7

#endif /* PARAMETER_HEADER_DD_H_ */
