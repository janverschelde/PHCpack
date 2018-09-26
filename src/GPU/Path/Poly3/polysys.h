// polysys.h contains prototypes of templated data types for polynomial
// systems with complex coefficients of different precisions.
// The corresponding definitions of the functions are in polysys.tpp.

#ifndef __POLYSYS_H__
#define __POLYSYS_H__

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

#include "polymon.h"
#include "polyeq.h"
#include "int_idx.h"

using namespace std;

// Polynomial system class
/*
   Read polynomial system from string or string arrays,
   and build array of polynomial equations, PolyEq.
   @sa PolyMon, PolyEq.
*/

template <class ComplexType, class RealType>
class PolySys
{
   public:

      int n_eq; // Number of Equations
      int dim;  // Dimension
      string* pos_var;
      vector<PolyEq<ComplexType,RealType>*> eq;
      // Array of polynomial equations
	
      PolyEq<ComplexType, RealType>* eq_space;

      int* max_deg_base;
      bool eval_base;

      int* job_number_level;
      int level;

      PolySys() // Constructor
      {
         n_eq = 0;
         dim = 0;
         pos_var = NULL;
         eq_space = NULL;
         max_deg_base = NULL;
         eval_base = false;
         job_number_level = NULL;
         level = 0;
      }

      ~PolySys() // Destructor
      {
         // delete eq_space;
         // delete pos_var;
         if(max_deg_base != NULL) delete[] max_deg_base;
         // cout << "sys destructed" << endl;
      }

      void read ( const string* sys_string, int n_eq, VarDict& pos_dict );
      // Reads a polynomial system from an array of strings.
      /*
        Read equation from string array, sys_string, using PolyEq::read().
        @param sys_string polynomial system string
        @param n number of equations or number of strings
        in string array sys_string
        @param pos_dict the dictionary of variables and their positions
        @sa PolyEq::read()
       */

      void read ( const string& sys_string, VarDict& pos_dict );
      // Reads a polynomial system from a string.
      /*
         First split the string by ';',
         then use PolyEq::read() to read each equation.
         @param sys_string polynomial system string
         @param pos_dict the dictionary of variables and their positions
         @sa PolyEq::read()
       */

      void read_file ( const string& file_name );
      // Reads a polynomial system from file.
      /*
         Reads the dimension on the first line and then the system.
         @param filename file name of polynomial system.
         File requirement: first line is dimension
         and the equations are separated by ';' symbols.
         @param pos_dict position dictionary for variables.
         If it is empty, the positions will be created by the order 
         of variables appearance first time in the file.
         @sa read_file(ifstream& myfile, Dict& pos_dict)
         Note that there should not be any blank lines
         between the lines that contain polynomials.
       */

      void read_file ( ifstream& myfile, VarDict& pos_dict );
      // Reads a polynomial system from file.
      /*
         Read dimension on the first line and then system.
         @param myfile file stream of polynomial system.
         File requirement: first line is dimension
         and the equations are separated by ';' symbols.
         @param pos_dict position dictionary for variables.
         If it is empty, the positions will be created by the order 
         of variables appearance first time in the file.
         @sa read_file(string file_name, Dict& pos_dict)
       */

      ComplexType* eval ( const ComplexType* x_val );
      // Evaluates a polynomial system at the point x_val.
      /*
         @param x_val ComplexType array of variables' values
         @return ComplexType array of system value
         On input, x_val points to as many complex numbers
         as the dimension of the system.
         On return is a pointer to a newly allocated array
         of as many complex numbers as the number of equations.
      */

      ComplexType* eval ( const ComplexType* x_val, ComplexType** deri );
      // Evaluates a polynomial system and its derivatives.
      /*
         Evaluate monomials and their derivatives by PolyMon::eval(),
         and add them up.
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return equation value
       */

      void eval ( const ComplexType* x_val, ComplexType* f_val,
                  ComplexType** deri_val );

      void print();
      // Prints the polynomial system.
      /*
         Use PolyEq::print() to print each equation in the system.
         @param pos_dict the dictionary of variables and their positions
         @return void
         @sa PolyEq::print() PolyMon::print()
       */
    
      void gpu_mon
         ( int& dim, int& level, int& workspace_size, int*& workspace_level,
           int*& n_mon_level, int& pos_size, unsigned short*& pos,
           int*& pos_level, int& sum_level, int*& n_sum_level,
           int& total_n_sum, int& sum_array_size, int*& sum_start,
           int*& sum_array );

      int workspace_size_block ( int start_level, int factor_size );

      int memory_size ( int factor_size );

      void job_number();

      void print_level();

      int job_number_block ( int start_level );

      void update_max_deg_base();

      ComplexType** eval_deg ( const ComplexType* x_val );

      void balance_eq ( const ComplexType* x_val );
};

#include "polysys.tpp"

#endif
