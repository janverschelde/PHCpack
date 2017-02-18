#include "parameter.h"

using namespace std;

void Parameter::write ( void )
{
   cout << " 1. Maximum number of steps                   : "
        << this->max_step << endl;
   cout << " 2. Number of points in the predictor         : "
        << this->n_predictor << endl;
   cout << " 3. Increase factor on the step size          : "
        << this->step_increase << endl;
   cout << " 4. Decrease factor on the step size          : "
        << this->step_decrease << endl;
   cout << " 5. Maximal step size along a path            : "
        << this->max_delta_t << endl;
   cout << " 6. Maximal step size at the end of a path    : "
        << this->max_delta_t_end << endl;
   cout << " 7. Minimum step size along a path            : "
        << this->min_delta_t << endl;
   cout << " 8. Tolerance on the residual                 : "
        << this->err_max_res << endl;
   cout << " 9. Tolerance on the corrector update         : "
        << this->err_max_delta_x << endl;
   cout << "10. Tolerance on the first corrector update   : "
        << this->err_max_first_delta_x << endl;
   cout << "11. Maximum number of Newton iterations       : "
        << this->max_it << endl;
   cout << "12. Tolerance for Newton's corrector method   : "
        << this->err_min_round_off << endl;
   cout << "13. Maximum number of Newton refinement steps : "
        << this->max_it_refine << endl;
   cout << "14. Tolerance for Newton's refinement method  : "
        << this->err_min_round_off_refine << endl;
}

void Parameter::set_value ( int idx, double val )
{
   if(idx == 1)
       this->max_step = (int) val;
   else if(idx == 2)
       this->n_predictor = (int) val;
   else if(idx == 3)
       this->step_increase = val;
   else if(idx == 4)
       this->step_decrease = val;
   else if(idx == 5)
       this->max_delta_t = val;
   else if(idx == 6)
       this->max_delta_t_end = val;
   else if(idx == 7)
       this->min_delta_t = val;
   else if(idx == 8)
       this->err_max_res = val;
   else if(idx == 9)
       this->err_max_delta_x = val;
   else if(idx == 10)
       this->err_max_first_delta_x = val;
   else if(idx == 11)
       this->max_it = (int) val;
   else if(idx == 12)
       this->err_min_round_off = val;
   else if(idx == 13)
       this->max_it_refine = (int) val;
   else if(idx == 14)
       this->err_min_round_off_refine = val;
}

void Parameter::get_value ( int idx, double* val )
{
   if(idx == 1)
       *val = (double) this->max_step;
   else if(idx == 2)
       *val = (double) this->n_predictor;
   else if(idx == 3)
       *val = this->step_increase;
   else if(idx == 4)
       *val = this->step_decrease;
   else if(idx == 5)
       *val = this->max_delta_t;
   else if(idx == 6)
       *val = this->max_delta_t_end;
   else if(idx == 7)
       *val = this->min_delta_t;
   else if(idx == 8)
       *val = this->err_max_res;
   else if(idx == 9)
       *val = this->err_max_delta_x;
   else if(idx == 10)
       *val = this->err_max_first_delta_x;
   else if(idx == 11)
       *val = (double) this->max_it;
   else if(idx == 12)
       *val = this->err_min_round_off;
   else if(idx == 13)
       *val = (double) this->max_it_refine;
   else if(idx == 14)
       *val = this->err_min_round_off_refine;
}

void Parameter::tune ( void )
{
   int idx = 1;

   while(idx > 0)
   {
      this->write();

      cout << "Type integer in 1..14 to change, 0 to exit : ";
      cin >> idx;

      if(idx > 0 && idx < 14)
      {
         cout << "Give a value for parameter " << idx << " : ";
         double val; cin >> val;
         this->set_value(idx,val);
      }
   }
}
