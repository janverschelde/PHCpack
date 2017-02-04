/* parameter.h, created on Dec 25, 2014 by yxc and edited by jv */

#ifndef PARAMETER_H_
#define PARAMETER_H_

// #include "Parameter_Header.h"
// the following definition is copied from Parameter_Header.h
#define MON_EVAL_METHOD 2
// 0 : reverse mode
// 1 : reverse mode with aligned memory for instructions
// 2 : tree mode
// 3 : for multiple evaluation (chosen when #paths > 1)

// Default values for parameters
#define N_PREDICTOR           4

#define MAX_STEP              400
#define MAX_DELTA_T           1E-1
#define MAX_DELTA_T_END       1E-2
#define MIN_DELTA_T           1E-7

#define MAX_IT                3
#define ERR_MIN_ROUND_OFF     1E-9

#define MAX_IT_REFINE                   5
#define ERR_MIN_ROUND_OFF_REFINE    1E-11

#define ERR_MAX_RES           1E-2
#define ERR_MAX_DELTA_X       1E-1
#define ERR_MAX_FIRST_DELTA_X 1E-2

#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7

#include <iostream>

class Parameter
{
   public:

      int n_predictor;

      int max_step;

      int max_it;
      double max_delta_t;
      double max_delta_t_end;
      double min_delta_t;
      double err_max_res;
      double err_max_delta_x;
      double err_max_first_delta_x;
      double err_min_round_off;
      int max_it_refine;
      double err_min_round_off_refine;
      double step_increase;
      double step_decrease;

      Parameter
       ( int n_predictor, int max_step, int max_it,
         double max_delta_t, double max_delta_t_end, double min_delta_t,
         double err_max_res, double err_max_delta_x,
         double err_max_first_delta_x, double err_min_round_off,
         double max_it_refine, double err_min_round_off_refine,
         double step_increase, double step_decrease )
      {
         this->n_predictor = n_predictor;
         this->max_step = max_step;
         this->max_it = max_it;
         this->max_delta_t = max_delta_t;
         this->max_delta_t_end = max_delta_t_end;
         this->min_delta_t = min_delta_t;
         this->err_max_res = err_max_res;
         this->err_max_delta_x = err_max_delta_x;
         this->err_max_first_delta_x = err_max_first_delta_x;
         this->err_min_round_off = err_min_round_off;
         this->max_it_refine = max_it_refine;
         this->err_min_round_off_refine = err_min_round_off_refine;
         this->step_increase = step_increase;
         this->step_decrease = step_decrease;
      }

      Parameter ( void ) // sets the default values for the parameters
      {
         this->n_predictor = N_PREDICTOR;
         this->max_step = MAX_STEP;
         this->max_it = MAX_IT;
         this->max_delta_t = MAX_DELTA_T;
         this->max_delta_t_end = MAX_DELTA_T_END;
         this->min_delta_t = MIN_DELTA_T;
         this->err_max_res = ERR_MAX_RES;
         this->err_max_delta_x = ERR_MAX_DELTA_X;
         this->err_max_first_delta_x = ERR_MAX_FIRST_DELTA_X;
         this->err_min_round_off = ERR_MIN_ROUND_OFF;
         this->max_it_refine = MAX_IT_REFINE;
         this->err_min_round_off_refine = ERR_MIN_ROUND_OFF_REFINE;
         this->step_increase = STEP_INCREASE;
         this->step_decrease = STEP_DECREASE;
      }
};

#endif /* PARAMETER_H_ */
