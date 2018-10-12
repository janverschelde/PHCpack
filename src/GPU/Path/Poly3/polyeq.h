// The file polyeq.h contains prototypes of templated data types for 
// polynomials with complex coefficients of different precisions.
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

/*
 * The polynomial equation class defines the representation and
 * the evaluation and differentiation methods for polynomials
 * with complex coefficients (of various precisions) in several variables.
 */

template <class ComplexType, class RealType>
class PolyEq
{
   public:

      ComplexType constant; // the constant of the polynomial
      int n_mon;            // number of monomials in the vector mon
      int dim;              // maximum number of variables in each monomial
      vector<PolyMon<ComplexType,RealType>*> mon;
      // monomial array of the equation of size n_mon
      // with nonzero constant, the polynomial has thus n_mon+1 terms

      int level;
      int* job_number_level;

      PolyEq() // returns the zero polynomial
      {
         constant.init(0.0,0.0);
         n_mon = 0;
         dim = 0;
         level = 0;
         job_number_level = NULL;
      }

      PolyEq ( int dim ) // returns the zero polynomial and sets dim
      {
         constant.init(0.0, 0.0);
         n_mon = 0;
         this -> dim = dim;
         level = 0;
         job_number_level = NULL;
      }

      ~PolyEq() // destructor
      {
          for(typename vector< PolyMon<ComplexType,RealType>* >::iterator
              it=mon.begin(); it<mon.end(); it++)
          {
             delete *it;
          }
      }

      void read ( const string& eq_string, VarDict& pos_dict, int verbose=0 );
      /*
       * Reads a polynomial from a string.
       * The dictionary pos_dict defines the positions of the variables.
       * If verbose, then extra output is written to screen.
       */

      void read ( const string& eq_string, VarDict& pos_dict,
                  int start, int end, int verbose=0 );
      /*
       * Reads a polynomial from a string, starting at position start
       * and ending at position end.
       * The dictionary pos_dict defines the positions of the variables.
       * If verbose, then extra output is written to screen.
       */

      void print ( const string* pos_var );
      /*
       * Writes the polynomial term by term, with the constant last.
       * The string pos_var defines the dictionary of the variables
       * and their positions.
       */

      ComplexType eval ( const ComplexType* x_val );
      /*
       * Evaluates the polynomial term by term,
       * using the straighforward algorithm to evaluate a monomial.
       * The input variable x_val must contain dim complex numbers.
       */
	
      ComplexType eval
       ( const ComplexType* x_val, ComplexType* deri, ComplexType *monderi );
      /*
       * Applies the algorithm of Speelpenning to evaluate and differentiate
       * all monomials in the polynomial at the values in x_val.
       * The function returns the function value and the values of the
       * derivatives are returned in the argument deri.
       * The dim values in monderi are auxiliary to store the intermediate
       * values of the derivatives of the monomials.
       * Results will be correct only if no variables appear with
       * power 2 or higher.
       */
	
      ComplexType eval ( const ComplexType* x_val, ComplexType* deri );
      /*
       * Applies the algorithm of Speelpenning to evaluate and differentiate
       * all monomials in the polynomial at the values in x_val.
       * The function returns the function value and the values of the
       * derivatives are returned in the argument deri.
       * Results will be correct only if no variables appear with
       * power 2 or higher.
       * This method allocates and deallocates the auxiliary data structure
       * monderi and wraps the previous eval method.
       * Because of this allocation and deallocation, this evaluation
       * method should not be used in a multithreaded application.
       */

      ComplexType eval
       ( const ComplexType* x_val, ComplexType* deri, ComplexType* monderi,
         ComplexType** deg_table );
      /*
       * Evaluates and differentiations the polynomial at a point
       * given by its coordinates in x_val and where all necessary 
       * powers of the coordinates are given in the parameter deg_table.
       * The value of the polynomial at the point is returned
       * and the values of all its partial derivatives at the point
       * are returned in the argument deri.
       * The argument monderi is an auxiliary array of the same size 
       * as x_val to hold all derivatives of the monomials.
       */

      ComplexType eval ( const ComplexType* x_val, ComplexType* deri,
                         ComplexType** deg_table );
      /*
       * Evaluates and differentiations the polynomial at a point
       * given by its coordinates in x_val and where all necessary 
       * powers of the coordinates are given in the parameter deg_table.
       * The value of the polynomial at the point is returned
       * and the values of all its partial derivatives at the point
       * are returned in the argument deri.
       * This method allocates and deallocates an auxiliary array
       * and should not be used in a multithreaded application.
       * Instead, use the previous eval with the monderi argument.
       */

      int memory_size ( int factor_size );

      void print_level();

      int workspace_size_block ( int start_level, int factor_size );

      int job_number_block ( int start_level );

      void update_max_deg ( int* max_deg );
      /*
       * Updates the array of integers in max_deg with the value
       * of the largest exponent, relative to its position.
       */
};

#include "polyeq.tpp"

#endif
