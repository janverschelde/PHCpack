// This is the main test on polynomial evaluation and differentiation.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "gqd_qd_util.h"

using namespace std;

template<class realD, class realH>
int ped_test 
 ( int prc, int BS, int dim, int NM, int NV, int deg, int r, int mode );
/*
 * Runs the test on polynomial evaluation and differentiation.
 *
 * ON ENTRY :
 *   prc      the working precision, 0 for d, 1 for dd, 2 for qd;
 *   BS       number of threads in a block, the block size;
 *   dim      dimension of the problem;
 *   NM       number of monomials;
 *   NV       number of variables;
 *   deg      highest degree of the variables;
 *   r        frequency of the runs;
 *   mode     mode of execution, 0, 1, or 2. */

int main ( int argc, char *argv[] )
{
   int prc,BS,dim,NM,NV,deg,r,mode;

   if(parse_args(argc,argv,&prc,&BS,&dim,&NM,&NV,&deg,&r,&mode) == 1) return 1;

   if(mode == 2)
      cout << "Testing polynomial evaluation and differentiation ...\n";

   if(prc == 0)
      ped_test<double, double>(prc,BS,dim,NM,NV,deg,r,mode);
   else if(prc == 1)
      ped_test<gdd_real, dd_real>(prc,BS,dim,NM,NV,deg,r,mode);
   else if(prc == 2)
      ped_test<gqd_real, qd_real>(prc,BS,dim,NM,NV,deg,r,mode);
   else
      cout << "Wrong level of precision " << prc
           << ", should be 0, 1, or 2." << endl;
}

template<class realD, class realH>
int ped_test 
 ( int prc, int BS, int dim, int NM, int NV, int deg, int r, int mode )
{
   if(mode == 2)
   {
      cout << "Generating a problem of dimension " << dim
           << " ..." << endl;
      cout << "number of monomials : " << NM << endl;
      cout << "number of variables : " << NV << endl;
      cout << "degree of each variable : " << deg << endl;
   }

   int pos_arr_h_int[NM*NV];
   int exp_arr_h_int[NM*NV];
   char pos_arr_h_char[NM*NV];
   char exp_arr_h_char[NM*NV];
   int ncoefs = NM*(NV+1);

   complexD<realD> *c_d = new complexD<realD>[ncoefs];
   complexH<realH> *c_h = new complexH<realH>[ncoefs];

   generate_system<realD, realH>
      (dim,NM,NV,deg,pos_arr_h_int,pos_arr_h_char,
       exp_arr_h_int,exp_arr_h_char,c_d,c_h);

   if(mode == 2) write_system<realH>(dim,NM,NV,c_h,pos_arr_h_int,exp_arr_h_int);
}
