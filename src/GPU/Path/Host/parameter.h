/* parameter.h, created on Dec 25, 2014 by yxc and edited by jv */

#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "Parameter_Header.h"

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
         // std::cout << "path_parameter.err_min_round_off = "
         //           << this->err_min_round_off << std::endl;
      }
};

#endif /* PARAMETER_H_ */
