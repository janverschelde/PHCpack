// Tests Newton's method with algorithmic differentiation.

#include <iostream>
#include <string>
#include "complexH.h"
#include "polysol.h"
#include "eval_host.h"
#include "path_host.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( string targetfile, string startfile );
// reads target and start system from file, calls path_test

template <class ComplexType, class RealType>
int path_test
 ( PolySys<ComplexType,RealType>& targetpols,
   PolySys<ComplexType,RealType>& startpols );
// Calls the path tracker on a homotopy defined by the target system
// in targetpols and the start system in startpols.

string* x_var ( string x, int dim );
// Generates a string of variables of size dim.

template <class ComplexType, class RealType>
bool read_homotopy_from_file
 ( PolySys<ComplexType,RealType>& targetsys,
   PolySys<ComplexType,RealType>& startsys,
   PolySolSet<ComplexType,RealType>* sols,
   string startname, string targetname );
// Reads a target system, start system, and solutions from file.

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
   cout << "Testing the path tracker ..." << endl;

   string name4target; // = "/Users/jan/PHCv2/Demo/cyclic5";
   cout << "\n-> give a file name for the target system : ";
   cin >> name4target;

   string name4start;
   cout << "\n-> give a file name for the start system : ";
   cin >> name4start;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   int fail = 0;

   if(choice == '0')
      fail = test<complexH<double>,double>(name4target,name4start);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name4target,name4start);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name4target,name4start);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int test ( string targetfile, string startfile )
{
   PolySys<ComplexType,RealType> targetsys;
   PolySys<ComplexType,RealType> startsys;
   // PolySolSet<ComplexType,RealType> sols;

   cout << "\nReading from file " << targetfile << " ..." << endl;
   targetsys.read_file(targetfile);
   cout << "The polynomials on the file " << targetfile << " :" << endl;
   targetsys.print();

   cout << "\nReading from file " << startfile << " ..." << endl;
   startsys.read_file(startfile);
   cout << "The polynomials on the file " << startfile << " :" << endl;
   startsys.print();

   ifstream solfile(startfile.c_str());
   PolySolSet<ComplexType,RealType> sols(solfile);
   cout << "The solutions on the file " << startfile << " :" << endl;
   sols.print();

   //bool success = read_homotopy_from_file<ComplexType,RealType>
   //   (targetsys,startsys,&sols,targetfile,startfile);

   return path_test<ComplexType,RealType>(targetsys,startsys);
}

string* x_var ( string x, int dim )
{
   string* var = new string[dim];

   for(int i=0; i< dim; i++)
   {
      ostringstream ss;
      ss << x << i;
      var[i] = ss.str();
   }
   return var;
}

template <class ComplexType, class RealType>
bool read_homotopy_from_file
 ( PolySys<ComplexType,RealType>& targetsys,
   PolySys<ComplexType,RealType>& startsys,
   PolySolSet<ComplexType,RealType>* sols,
   string startname, string targetname )
{
   ifstream targetfile(targetname.c_str());
   if(targetfile.is_open() == false)
   {
      cout << "Unable to open file with name \""
           << targetname << "\"." << endl;
      return false;
   }
   ifstream startfile(startname.c_str());
   if(startfile.is_open() == false)
   {
      cout << "Unable to open file with name \""
           << startname << "\"." << endl;
      return false;
   }
   VarDict pos_dict;

   startsys.read_file(startfile,pos_dict);

   int dim = startsys.dim;

   string x_name = "x";
   string* x_names = x_var(x_name, dim);
   startsys.pos_var = x_names;

   cout << "the start system read from file" << endl;
   startsys.print();

   // if(sols != NULL) sols -> init(targetfile, pos_dict);

   targetsys.read_file(targetfile,pos_dict);

   string x_name_target = "x";
   string* x_names_target = x_var(x_name_target,dim);
   targetsys.pos_var = x_names_target;

   cout << "the target system read from file" << endl;
   targetsys.print();

   return true;
}

template <class ComplexType, class RealType>
int path_test
 ( PolySys<ComplexType,RealType>& targetpols,
   PolySys<ComplexType,RealType>& startpols )
{
   cout << "The dimension : " << targetpols.dim << endl;

   CPUInstHom<ComplexType,RealType> ped; // data for eval_host
   Workspace<ComplexType> wrk;
   ComplexType alpha,t;
   alpha.init(0.0,0.0); // initialize the data for eval_host
   t.init(0.0,0.0);
   ped.init(targetpols,startpols,targetpols.dim,targetpols.n_eq,0,alpha);
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
