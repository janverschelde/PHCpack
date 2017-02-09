// Test on the interactive tuning of parameters.

#include "parameter.h"

using namespace std;

void write ( Parameter &pars );
// Writes the values in pars to screen.

void set_value ( Parameter &pars, int idx, double val );
// Sets the value of the parameter with idx to the value val.
// Required: idx is in the range from 1 to 14.

void tune ( Parameter &pars );
// Interactive tuning of the parameters in pars.

int main ( void )
{
   int precision;
  
   cout << "Give the precision (16, 32, or 64) : ";
   cin >> precision;

   Parameter pars(precision);

   tune(pars);

   return 0;
}

void tune ( Parameter &pars )
{
   int idx = 1;

   while(idx > 0)
   {
      write(pars);

      cout << "Type integer in 1..14 to change, 0 to exit : ";
      cin >> idx;

      if(idx > 0 && idx < 14)
      {
         cout << "Give a value for parameter " << idx << " : ";
         double val; cin >> val;
         set_value(pars,idx,val);
      }
   }
}

void write ( Parameter &pars )
{
   cout << " 1. Maximum number of steps                   : "
        << pars.max_step << endl;
   cout << " 2. Number of points in the predictor         : "
        << pars.n_predictor << endl;
   cout << " 3. Increase factor on the step size          : "
        << pars.step_increase << endl;
   cout << " 4. Decrease factor on the step size          : "
        << pars.step_decrease << endl;
   cout << " 5. Maximal step size along a path            : "
        << pars.max_delta_t << endl;
   cout << " 6. Maximal step size at the end of a path    : "
        << pars.max_delta_t_end << endl;
   cout << " 7. Minimum step size along a path            : "
        << pars.min_delta_t << endl;
   cout << " 8. Tolerance on the residual                 : "
        << pars.err_max_res << endl;
   cout << " 9. Tolerance on the corrector update         : "
        << pars.err_max_delta_x << endl;
   cout << "10. Tolerance on the first corrector update   : "
        << pars.err_max_first_delta_x << endl;
   cout << "11. Maximum number of Newton iterations       : "
        << pars.max_it << endl;
   cout << "12. Tolerance for Newton's corrector method   : "
        << pars.err_min_round_off << endl;
   cout << "13. Maximum number of Newton refinement steps : "
        << pars.max_it_refine << endl;
   cout << "14. Tolerance for Newton's refinement method  : "
        << pars.err_min_round_off_refine << endl;
}

void set_value ( Parameter &pars, int idx, double val )
{
   if(idx == 1)
       pars.max_step = (int) val;
   else if(idx == 2)
       pars.n_predictor = (int) val;
   else if(idx == 3)
       pars.step_increase = val;
   else if(idx == 4)
       pars.step_decrease = val;
   else if(idx == 5)
       pars.max_delta_t = val;
   else if(idx == 6)
       pars.max_delta_t_end = val;
   else if(idx == 7)
       pars.min_delta_t = val;
   else if(idx == 8)
       pars.err_max_res = val;
   else if(idx == 9)
       pars.err_max_delta_x = val;
   else if(idx == 10)
       pars.err_max_first_delta_x = val;
   else if(idx == 11)
       pars.max_it = (int) val;
   else if(idx == 12)
       pars.err_min_round_off = val;
   else if(idx == 13)
       pars.max_it_refine = (int) val;
   else if(idx == 14)
       pars.err_min_round_off_refine = val;
}
