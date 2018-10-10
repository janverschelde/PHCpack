// Tests Newton's method with algorithmic differentiation.

#include <iostream>
#include <iomanip>
#include <string>
#include "complexH.h"
#include "eval_host.h"
#include "polysolset.h"
#include "newton_host.h"

using namespace std;

template <class ComplexType, class RealType>
int newton_test
 ( PolySys<ComplexType,RealType>& polynomials,
   PolySolSet<ComplexType,RealType>& solutions, int precision );
/*
 * Calls Newton's method on the polynomials, starting at the values in
 * solutions, in double, double double, or quad double precision, 
 * for precision respectively equal to 16, 32, or 64. */

template <class ComplexType, class RealType>
int test ( string filename, int precision );
/*
 * Reads a polynomial system and solutions from file, calls newton_test,
 * where the value of precision should be 16, 32, or 62,
 * respectively for double, double double, or quad double precision. */

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
      fail = test<complexH<double>,double>(name,16);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name,32);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name,64);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int test ( string filename, int precision )
{
   PolySys<ComplexType,RealType> polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   ifstream solfile(filename.c_str());
   PolySolSet<ComplexType,RealType> sols(solfile);
   cout << "The solutions on the file " << filename << " :" << endl;
   sols.print();

   return newton_test<ComplexType,RealType>(polynomials,sols,precision);
}

template <class ComplexType, class RealType>
int newton_test
 ( PolySys<ComplexType,RealType>& polynomials,
   PolySolSet<ComplexType,RealType>& solutions, int precision )
{
   cout << "The dimension : " << polynomials.dim << endl;

   ComplexType* arg = new ComplexType[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(i,0.0);

   CPUInstHom<ComplexType,RealType> ped; // data for eval_host
   Workspace<ComplexType> wrk;
   ComplexType alpha,t;
   alpha.init(0.0,0.0); // initialize the data for eval_host
   t.init(0.0,0.0);
   ped.init(polynomials,polynomials.dim,polynomials.n_eq,0,alpha,1);

   cout << "Number of equations n_eq in ped : " << ped.n_eq << endl;

   ped.init_workspace(wrk);
   Parameter pars(precision);

   pars.tune();

   double tSec_Eval = 0;
   double tSec_MGS = 0;

   wrk.sol = solutions.get_sol(0);
   wrk.update_x_t(wrk.sol,t);

   cout << "The coordinates of the first solution :" << endl;
   for(int i=0; i<polynomials.dim; i++)
      cout << setprecision(precision) << wrk.sol[i];

   bool success_cpu
      //= CPU_Newton_Refine<ComplexType,RealType>
      = CPU_Newton<ComplexType,RealType>
           (wrk,ped,pars,tSec_Eval,tSec_MGS);

   cout << "The coordinates after Newton's method :" << endl;
   for(int i=0; i<polynomials.dim; i++)
      cout << setprecision(precision) << wrk.x[i];

   return 0;
}
