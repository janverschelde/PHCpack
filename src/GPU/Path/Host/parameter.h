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
#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7
#define MAX_DELTA_T           1E-1
#define MAX_DELTA_T_END       1E-2
#define MIN_DELTA_T           1E-7
#define ERR_MAX_RES           1E-6
#define ERR_MAX_DELTA_X       1E-6
#define ERR_MAX_FIRST_DELTA_X 1E-2

#define D_MAX_STEP            1000
#define DD_MAX_STEP           2000
#define QD_MAX_STEP           3000

#define D_MAX_IT              3
#define DD_MAX_IT             4
#define QD_MAX_IT             5
#define D_ERR_MIN_ROUND_OFF   1E-9
#define DD_ERR_MIN_ROUND_OFF  1E-14
#define QD_ERR_MIN_ROUND_OFF  1E-26

#define D_MAX_IT_REFINE                  3
#define DD_MAX_IT_REFINE                 4
#define QD_MAX_IT_REFINE                 5
#define D_ERR_MIN_ROUND_OFF_REFINE   1E-11
#define DD_ERR_MIN_ROUND_OFF_REFINE  1E-22
#define QD_ERR_MIN_ROUND_OFF_REFINE  1E-40

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

      Parameter ( int prc ) // sets the defaults for precision prc
      {
         this->n_predictor = N_PREDICTOR;
         this->step_increase = STEP_INCREASE;
         this->step_decrease = STEP_DECREASE;
         this->max_delta_t = MAX_DELTA_T;
         this->max_delta_t_end = MAX_DELTA_T_END;
         this->min_delta_t = MIN_DELTA_T;
         this->err_max_res = ERR_MAX_RES;
         this->err_max_delta_x = ERR_MAX_DELTA_X;
         this->err_max_first_delta_x = ERR_MAX_FIRST_DELTA_X;
         if(prc <= 16)
         {
            this->max_step = D_MAX_STEP;
            this->max_it = D_MAX_IT;
            this->err_min_round_off = D_ERR_MIN_ROUND_OFF;
            this->max_it_refine = D_MAX_IT_REFINE;
            this->err_min_round_off_refine = D_ERR_MIN_ROUND_OFF_REFINE;
         }
         else if(prc <= 32)
         {
            this->max_step = DD_MAX_STEP;
            this->max_it = DD_MAX_IT;
            this->err_min_round_off = DD_ERR_MIN_ROUND_OFF;
            this->max_it_refine = DD_MAX_IT_REFINE;
            this->err_min_round_off_refine = DD_ERR_MIN_ROUND_OFF_REFINE;
         }
         else
         {
            this->max_step = QD_MAX_STEP;
            this->max_it = QD_MAX_IT;
            this->err_min_round_off = QD_ERR_MIN_ROUND_OFF;
            this->max_it_refine = QD_MAX_IT_REFINE;
            this->err_min_round_off_refine = QD_ERR_MIN_ROUND_OFF_REFINE;
         }
      }
     
      void write ( void );
      // Writes the values in pars to screen.

      void set_value ( int idx, double val );
      // Sets the value of the parameter with idx to the value val.
      // Required: idx is in the range from 1 to 14.

      void get_value ( int idx, double* val );
      // Gets the value of the parameter with idx 
      // and assigns the value to the parameter val.
      // Required: idx is in the range from 1 to 14.

      void tune ( void );
      // Interactive tuning of the parameters in pars.
};

#endif /* PARAMETER_H_ */
