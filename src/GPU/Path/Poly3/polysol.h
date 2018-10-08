// Defines solutions of polynomial systems over complex numbers.

#ifndef __POLYSOL_H__
#define __POLYSOL_H__

#include <iostream>
#include <fstream>
#include "stdlib.h"
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include "linklist.h"
#include "dict.h"
#include "utilities.h"

template <class ComplexType, class RealType>
class PolySol
{
   public:

      int dim;          // solution number
      ComplexType* sol; // coordinates of the solution point

      int idx;
      int path_idx;
      int m;  // multiplicity

      ComplexType t; // value of the continuation parameter
      RealType err;  // forward error
      RealType rco;  // inverse of condition number estimate
      RealType res;  // residual, backward error

      string info; // Success / Fail / Infinity

      PolySol()
      {
         dim = 0;
         idx = 0;
         path_idx = 0;
         m = 0;
         t = ComplexType(0.0,0.0);
         sol = NULL;
         err = 0.0;
         rco = 0.0;
         res = 0.0;
      }

      void init ( ifstream& myfile, int dim );

      void init ( ifstream& myfile, int dim, VarDict& pos_dict );

      void init
       ( ComplexType* sol, int dim, RealType max_residual, 
         RealType max_delta_x, int path_idx, string path_info );

      void init ( int dim, RealType t_real, RealType t_imag, RealType* sol,
                  RealType max_delta_x, RealType rco, RealType max_residual,
                  int m, int path_idx, string path_info );

      PolySol ( ifstream& myfile, int dim )
      {
         init(myfile, dim);
      }

      PolySol ( ifstream& myfile, int dim, VarDict& pos_dict )
      {
         init(myfile, dim, pos_dict);
      }

      PolySol ( ComplexType* sol, int dim, RealType max_residual = 0,
                RealType max_delta_x=0,
                int path_idx=0, string path_info="" )
      {
         init(sol,dim,max_residual,max_delta_x,path_idx,path_info);
      }

      PolySol ( int dim, RealType t_real, RealType t_imag, RealType* sol, 
                RealType max_delta_x=0, RealType rco=0,
                RealType max_residual = 0,
                int m=0, int path_idx=0, string path_info="" )
      { 
         init(dim,t_real,t_imag,sol,max_delta_x,rco,max_residual,m,
              path_idx,path_info);
      }

      ~PolySol()
      {
         delete[] sol;
      }

      bool operator == ( const PolySol<ComplexType,RealType>& that );

      bool operator< ( PolySol<ComplexType,RealType>& that );

      void print();

      void print_short();

      void print_info();

      void print_info ( string* pos_var );

      ComplexType* get_sol();
};

template <class ComplexType, class RealType>
bool compare_sol 
 ( PolySol<ComplexType,RealType>* sol1, PolySol<ComplexType,RealType>* sol2 );
// compares the solutions sol1 with sol2

#include "polysol.tpp"

#endif
