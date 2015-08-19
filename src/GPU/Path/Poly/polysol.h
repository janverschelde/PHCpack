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
//#include "job.h"
#include "dict.h"
#include "utilities.h"

class PolySol
{
   public:

      int dim; // solution number
      CT* sol;

      int idx;
      int path_idx;
      int m;  // multiplicity

      CT t;
      T1 err; // forward error
      T1 rco; // inverse of condition number
      T1 res; // residual, backward error

      string info; // Success / Fail / Infinity

      PolySol()
      {
         dim = 0;
         idx = 0;
         path_idx = 0;
         m = 0;
         t = CT(0.0,0.0);
         sol = NULL;
         err = 0.0;
         rco = 0.0;
         res = 0.0;
      }

      void init ( ifstream& myfile, int dim );

      void init ( ifstream& myfile, int dim, VarDict& pos_dict );

      void init
       ( CT* sol, int dim, T1 max_residual, T1 max_delta_x, int path_idx,
         string path_info );

      void init ( int dim, T1 t_real, T1 t_imag, T1* sol,
                  T1 max_delta_x, T1 rco, T1 max_residual,
                  int m, int path_idx, string path_info );

      PolySol ( ifstream& myfile, int dim )
      {
         init(myfile, dim);
      }

      PolySol ( ifstream& myfile, int dim, VarDict& pos_dict )
      {
         init(myfile, dim, pos_dict);
      }

      PolySol ( CT* sol, int dim, T1 max_residual = 0, T1 max_delta_x=0,
                int path_idx=0, string path_info="" )
      {
         init(sol,dim,max_residual,max_delta_x,path_idx,path_info);
      }

      PolySol ( int dim, T1 t_real, T1 t_imag, T1* sol, 
                T1 max_delta_x=0, T1 rco=0, T1 max_residual = 0,
                int m=0, int path_idx=0, string path_info="" )
      { 
         init(dim,t_real,t_imag,sol,max_delta_x,rco,max_residual,m,
              path_idx,path_info);
      }

      ~PolySol()
      {
         delete[] sol;
      }

      bool operator == (const PolySol& that);

      bool operator<(PolySol& that);

      void print();

      void print_short();

      void print_info();

      void print_info(string* pos_var);

      CT* get_sol();
};

bool compare_sol(PolySol* sol1, PolySol* sol2);

class PolySolSet
{
   public:

      int n_sol;
      int dim;
      vector<PolySol*> sols;

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

      bool find_same_sol ( PolySol* tmp_sol );

      int count_same_sol ( PolySol* tmp_sol );

      void add_sol
       ( CT* new_sol, T1 max_residual=0, T1 max_delta_x=0, int path_idx=0,
         string path_info="" );

      void add_sol ( PolySol* tmp_sol );

      void change_sol ( int idx, CT* coords );
        // updates the coordinates of the solution with index idx,
        // using the values in coords

      bool add_diff_sol ( CT* new_sol );

      void print();

      void print_info ( string* pos_var );

      void print_short();

      CT* get_sol(int idx);

      void sort_set();

      void compare ( PolySolSet& that );
};

#endif
