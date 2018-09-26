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

template <class ComplexType, class RealType>
class PolySolSet
{
   public:

      int n_sol;
      int dim;
      vector<PolySol<ComplexType,RealType>*> sols;

      void init ( ifstream& myfile );

      void init ( ifstream& myfile, VarDict& pos_dict );

      PolySolSet()
      {
         n_sol = 0;
         dim = 0;
      }

      PolySolSet ( ifstream& myfile )
      {
         init(myfile);
      }

      PolySolSet ( int dim )
      {
         this->dim = dim;
         n_sol = 0;
      }

      ~PolySolSet()
      {
         // std::cout << "Delete PolySolSet" << std::endl;
         for(int i=0; i<n_sol; i++)
         {
            delete sols[i];
         }
      }

      bool find_same_sol ( PolySol<ComplexType,RealType>* tmp_sol );

      int count_same_sol ( PolySol<ComplexType,RealType>* tmp_sol );

      void add_sol
       ( ComplexType* new_sol, RealType max_residual=0, 
         RealType max_delta_x=0, int path_idx=0, string path_info="" );

      void add_sol ( PolySol<ComplexType,RealType>* tmp_sol );

      void change_sol ( int idx, ComplexType* coords );
      // updates the coordinates of the solution with index idx,
      // using the values in coords

      void change_solt ( int idx, ComplexType* coords, ComplexType* tval );
      // updates the coordinates of the solution with index idx,
      // using the values in coords; also updates the t with tval

      bool add_diff_sol ( ComplexType* new_sol );

      void print();

      void print_info ( string* pos_var );

      void print_short();

      ComplexType* get_sol ( int idx );

      void get_solt ( int idx, ComplexType* sol, ComplexType* t );
      // returns in sol the coordinates of the solution with index idx
      // and in t its corresponding t value

      void sort_set();

      void compare ( PolySolSet<ComplexType,RealType>& that );
};

#include "polysol.tpp"

#endif
