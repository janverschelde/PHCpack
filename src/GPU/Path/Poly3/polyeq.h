// polyeq.h contains prototypes of templated data types for polynomial
// equations with complex coefficients of different precisions.
// The corresponding definitions of the functions are in polyeq.tpp.

#ifndef __POLYEQ_H__
#define __POLYEQ_H__

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

using namespace std;

// Polynomial equation class
/*
   Read polynomial equation from string,
   and build array of monomials, PolyMon.
   PolyEq is used to construct polynomial system, PolySys.
   @sa PolyMon PolySys
*/

template <class ComplexType, class RealType>
class PolyEq
{
   public:

      ComplexType constant; // Constant of equation
      int n_mon;   // Number of monomials
      int dim;     // Number of monomials
      vector<PolyMon<ComplexType,RealType>*> mon;
      // Monomial array of the equation

      int level;
      int* job_number_level;

      PolyEq() // Constructor 
      {
         constant.init(0.0,0.0);
         n_mon = 0;
         dim = 0;
         level = 0;
         job_number_level = NULL;
      }

      PolyEq(int dim) // Constructor 
      {
         constant.init(0.0, 0.0);
         n_mon = 0;
         this -> dim = dim;
         level = 0;
         job_number_level = NULL;
      }

      ~PolyEq() // Destructor
      {
          for(typename vector< PolyMon<ComplexType,RealType>* >::iterator
              it=mon.begin(); it<mon.end(); it++)
          {
             delete *it;
          }
      }

      void read ( const string& eq_string, VarDict& pos_dict );
      // Reads an equation from a string
      /*
         First, indentify monomials position.
         Then, construct an array of monomials, use PolyMon::read().
         @param eq_string equation string
         @param pos_dict the dictionary of variables and their positions
         @sa PolyMon::read()
       */

      void read ( const string& eq_string, VarDict& pos_dict,
                  int start, int end );
      // Reads an equation from certain part of a string
      /*
         First, indentify monomials position.
         Then, construct an array of monomials, use PolyMon::read().
         @param eq_string equation string
         @param pos_dict the dictionary of variables and their positions
         @param start the start position of equation in the string 
         @param end the end position of equation in the string
         @sa PolyMon::read()
       */

      void print ( const string* pos_var );
      // Print equation
      /*
         Print monomial by monomial in the equation.
         Constant comes the last.
         @param pos_dict the dictionary of variables and their positions
       */

      ComplexType eval ( const ComplexType* x_val );
      // Evaluate equation
      /*
         @param x_val ComplexType array of variables' values
         @return ComplexType equation value
       */
	
      ComplexType eval ( const ComplexType* x_val, ComplexType* deri );
      // Evaluate equation and its derivatives
      /*
         Evaluate monomials and their derivative by PolyMon::eval(),
         and add them up.
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return equation value
       */

      ComplexType eval ( const ComplexType* x_val, ComplexType* deri,
                         ComplexType** deg_table );

      int memory_size ( int factor_size );

      void print_level();

      int workspace_size_block ( int start_level, int factor_size );

      int job_number_block ( int start_level );

      void update_max_deg ( int* max_deg );
};

#include "polyeq.tpp"

#endif
