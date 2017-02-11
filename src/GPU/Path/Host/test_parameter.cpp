// Test on the interactive tuning of parameters.

#include <fstream>
#include "parameter.h"

using namespace std;

void write_to_file ( ofstream &outfile, Parameter &pars );
// Writes the parameters in pars to the file outfile,
// which must be opened for output.

int main ( void )
{
   int precision;
  
   cout << "Give the precision (16, 32, or 64) : ";
   cin >> precision;

   Parameter pars(precision);

   pars.tune();

   string filename;

   cout << "Give a file name : ";
   cin >> filename;

   ofstream outfile(filename.c_str());

   cout << "Writing to " << filename << " ..." << endl;;

   write_to_file(outfile,pars);

   outfile.close();

   return 0;
}

void write_to_file ( ofstream &outfile, Parameter &pars )
{
   outfile << " 1. Maximum number of steps                   : "
           << pars.max_step << endl;
   outfile << " 2. Number of points in the predictor         : "
           << pars.n_predictor << endl;
   outfile << " 3. Increase factor on the step size          : "
           << pars.step_increase << endl;
   outfile << " 4. Decrease factor on the step size          : "
           << pars.step_decrease << endl;
   outfile << " 5. Maximal step size along a path            : "
           << pars.max_delta_t << endl;
   outfile << " 6. Maximal step size at the end of a path    : "
           << pars.max_delta_t_end << endl;
   outfile << " 7. Minimum step size along a path            : "
           << pars.min_delta_t << endl;
   outfile << " 8. Tolerance on the residual                 : "
           << pars.err_max_res << endl;
   outfile << " 9. Tolerance on the corrector update         : "
           << pars.err_max_delta_x << endl;
   outfile << "10. Tolerance on the first corrector update   : "
           << pars.err_max_first_delta_x << endl;
   outfile << "11. Maximum number of Newton iterations       : "
           << pars.max_it << endl;
   outfile << "12. Tolerance for Newton's corrector method   : "
           << pars.err_min_round_off << endl;
   outfile << "13. Maximum number of Newton refinement steps : "
           << pars.max_it_refine << endl;
   outfile << "14. Tolerance for Newton's refinement method  : "
           << pars.err_min_round_off_refine << endl;
}
