// Tests the start system by evaluation of all its solutions.

#include <iostream>
#include <iomanip>
#include <string>
#include "complexH.h"
#include "polysys.h"
#include "polysolset.h"
#include "jobqueue.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( string startfile, int precision );
/*
 * Reads start system from file, calls start_test.
 * The precision should be 16, 32, or 64 for printing solutions.
 */

template <class ComplexType, class RealType>
int start_test
 ( PolySys<ComplexType,RealType>& startpols,
   PolySolSet<ComplexType,RealType>& sols, int precision );
/*
 * Evaluates all start solutions in sols in startpols.
 * The precision should be 16 (double), 32 (double double),
 * or 64 (quad double), for use in the printing of the solutions.
 */

template <class ComplexType, class RealType>
int crew_start_test
 ( int crewsize, PolySys<ComplexType,RealType>& startpols,
   PolySolSet<ComplexType,RealType>& sols, int precision );
/*
 * Evaluates all start solutions in sols in startpols.
 * The precision should be 16 (double), 32 (double double),
 * or 64 (quad double), for use in the printing of the solutions,
 * using a crew of crewsize threads.
 */

int main ( void )
{
   cout << "Testing the start system ..." << endl;

   string name4start;
   cout << "\n-> give a file name for the start system : ";
   cin >> name4start;

   char choice;

   cout << endl << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   int fail = 0;

   if(choice == '0')
      fail = test<complexH<double>,double>(name4start,16);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name4start,32);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name4start,64);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int start_test
 ( PolySys<ComplexType,RealType>& startpols,
   PolySolSet<ComplexType,RealType>& sols, int precision )
{
   const int neq = startpols.n_eq;
   const int dim = startpols.dim;

   ComplexType* funval = new ComplexType[neq];
   ComplexType** jacmat = new ComplexType*[neq];
   for(int idx=0; idx<neq; idx++) jacmat[idx] = new ComplexType[dim];

   ComplexType sum(0.0,0.0);

   for(int cnt=0; cnt<sols.n_sol; cnt++)
   {
      ComplexType* point = sols.get_sol(0);

      startpols.eval(point, funval, jacmat);
      for(int idx=0; idx<dim; idx++) sum = sum + funval[idx];
   }

   cout << "The sum of all function values : " << sum << endl;

   return 0;
}

template <class ComplexType, class RealType>
int crew_start_test
 ( int crewsize, PolySys<ComplexType,RealType>& startpols,
   PolySolSet<ComplexType,RealType>& sols, int precision )
{
   const int nbsols = sols.n_sol;
  
   ComplexType** points; // stores solution coordinates
   points = (ComplexType**) calloc(nbsols,sizeof(ComplexType*));

   for(int idx=0; idx<nbsols; idx++) points[idx] = sols.get_sol(idx);

   cout << scientific << setprecision(precision);

   for(int idx=0; idx<nbsols; idx++)
   {
      cout << "The coordinates of start solution " << idx
           << " :" << endl;
      for(int k=0; k<startpols.dim; k++) cout << points[idx][k];
   }

   JobQueue jobs(nbsols);
   jobs.write();

   return 0;
}

template <class ComplexType, class RealType>
int test ( string startfile, int precision )
{
   PolySys<ComplexType,RealType> startsys;

   cout << endl << "Reading from file " << startfile << " ..." << endl;
   startsys.read_file(startfile);
   cout << "The polynomials on the file " << startfile << " :" << endl;
   startsys.print();

   ifstream solfile(startfile.c_str());
   PolySolSet<ComplexType,RealType> sols(solfile);

   cout << endl << "-> read " << sols.n_sol
        << " solutions from " << startfile << endl;

   // cout << "The solutions on the file " << startfile << " :" << endl;
   // sols.print();

   cout << endl << "Give the number of threads (0 for no multithreading) : ";
   int nbthreads; cin >> nbthreads;

   int fail = 0;

   if(nbthreads == 0)
      fail = start_test<ComplexType,RealType>(startsys,sols,precision);
   else
      fail = crew_start_test<ComplexType,RealType>
                (nbthreads,startsys,sols,precision);

   return fail;
}
