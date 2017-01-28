// Tests Newton's method with algorithmic differentiation.

#include <iostream>
#include <string>
#include "complexH.h"
#include "eval_host.h"
#include "path_host.h"

using namespace std;

template <class ComplexType, class RealType>
int path_test ( PolySys<ComplexType,RealType>& polynomials );
// calls Newton's method on the polynomials

template <class ComplexType, class RealType>
int test ( string filename );
// reads a double precision system from file, calls path_test

// Parameters
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

int main ( void )
{
   cout << "Testing the polynomials ..." << endl;

   string name; // = "/Users/jan/PHCv2/Demo/cyclic5";
   cout << "-> give a file name : "; cin >> name;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   cout << "\nReading from file " << name << " ..." << endl;

   int fail = 0;

   if(choice == '0')
      fail = test<complexH<double>,double>(name);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int test ( string filename )
{
   PolySys<ComplexType,RealType> polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   return path_test<ComplexType,RealType>(polynomials);
}

template <class ComplexType, class RealType>
int path_test ( PolySys<ComplexType,RealType>& polynomials )
{
   cout << "The dimension : " << polynomials.dim << endl;

   ComplexType* arg = new ComplexType[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(i,0.0);

   CPUInstHom<ComplexType,RealType> ped; // data for eval_host
   Workspace<ComplexType> wrk;
   ComplexType alpha,t;
   alpha.init(0.0,0.0); // initialize the data for eval_host
   t.init(0.0,0.0);
   ped.init(polynomials,polynomials.dim,polynomials.n_eq,0,alpha);
   ped.init_workspace(wrk);

   Parameter pars(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T,
      MAX_DELTA_T_END, MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X,
      ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, MAX_IT_REFINE,
      ERR_MIN_ROUND_OFF_REFINE, STEP_INCREASE, STEP_DECREASE);

   double tSecPred = 0;
   double tSecEval = 0;
   double tSecMGS = 0;

   bool success_cpu = path_tracker<ComplexType,RealType>
      (wrk,ped,pars,tSecPred,tSecEval,tSecMGS);

   return 0;
}
