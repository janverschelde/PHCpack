// Defines sets of solutions of polynomial systems over complex numbers.

#ifndef __POLYSOLSET_H__
#define __POLYSOLSET_H__

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
#include "polysol.h"

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

#include "polysolset.tpp"

#endif
