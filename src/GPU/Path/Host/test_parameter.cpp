// Test on the interactive tuning of parameters.

#include "parameter.h"

using namespace std;

void write ( Parameter &pars );
// Writes the values in pars to screen.

int main ( void )
{
   int precision;
  
   cout << "Give the precision (16, 32, or 64) : ";
   cin >> precision;

   Parameter pars(precision);

   write(pars);

   return 0;
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
