// The file polysys.h contains prototypes of templated data types for
// polynomial systems with complex coefficients of different precisions.
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

/*
 * The class PolySys defines the representation of a polynomial system,
 * provides methods to read a symbolic formulation of the polynomials,
 * and methods to efficiently evaluate and differentiate the system
 * at a given point.
 * The coefficient type is complex and the templated construction
 * allows for different precisions.
 */

template <class ComplexType, class RealType>
class PolySys
{
   public:

      int n_eq; // the number of equations in the system
      int dim;  // the dimension equals the number of variables
      string* pos_var;
      vector<PolyEq<ComplexType,RealType>*> eq;
      // an array with the polynomial equations in the system

      PolyEq<ComplexType, RealType>* eq_space;

      int* max_deg_base; // computed by update_max_deg_base()
      bool eval_base;    // flag if there is a common factor to evaluate

      int* job_number_level;
      int level;

      PolySys() // constructor initializes to zero and to NULL
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

      ~PolySys() // the destructor deallocates occupied space
      {
         // delete eq_space;
         // delete pos_var;
         if(max_deg_base != NULL) delete[] max_deg_base;
         // cout << "sys destructed" << endl;
      }

      void read ( const string* sys_string, int n_eq, VarDict& pos_dict,
                  int verbose=0 );
      /*
       * Reads a polynomial systems from sys_string, an array of strings.
       * The number of strings equals n_eq, which equals the number
       * of polynomials in the system.
       * The method read() defined in the class PolyEq is applied.
       * The dictionary pos_dict defines the position of the variables.
       * If verbose > 0, then extra output is written to screen.
       */

      void read ( const string& sys_string, VarDict& pos_dict,
                  int verbose=0 );
      /*
       * Reads a polynomial system from the string sys_string.
       * Each polynomial in the string ends with a semicolon ';'.
       * The method read() defined in the class PolyEq is applied.
       * The dictionary pos_dict defines the positions of the variables.
       * If verbose > 0, then extra output is written to screen.
       */
 
      void read_file ( const string& file_name, int verbose=0 );
      /*
       * Reads the dimension of the system on the first line on file
       * and then reads the polynomial in the system.
       * The polynomials must end with a semicolon ';'.
       * There should not be any blank lines between the lines
       * that contain the polynomials.
       * If verbose > 0, then extra output is written to screen.
       */

      void read_file ( ifstream& myfile, VarDict& pos_dict, int verbose=0 );
      /*
       * Reads the dimension of the system on the first line on file
       * and then reads the polynomial in the system.
       * The polynomials must end with a semicolon ';'.
       * There should not be any blank lines between the lines
       * that contain the polynomials.
       * The dictionary pos_dict defines the position of the variables.
       * If the dictionary is empty, then the positions of the
       * variables are defined in the order of their appearance on file.
       * If verbose > 0, then extra output is written to screen.
       */

      ComplexType** allocate_deg_table ( void );
      /*
       * Allocates memory for the powers of the coordinates of the point
       * for evaluation of polynomials with higher degrees.
       * The allocated table is returned.
       */

      void compute_deg_table
       ( const ComplexType* x_val, ComplexType** deg_table );
      /*
       * Given in x_val the coordinates of a point and
       * an allocated memory space for the degree table in deg_table,
       * computes the powers of the variables at the coordinates in x_val.
       */

      ComplexType* eval ( const ComplexType* x_val );
      /*
       * Applies the straightforward algorithm to evaluate
       * the polynomial system at the values of the coordinates in x_val.
       * On input, x_val points to as many complex numbers
       * as the dimension of the system.
       * On return is a pointer to a newly allocated array
       * of as many complex numbers as the number of equations.
       * This function should only be used for testing purposes.
       */

      ComplexType* eval ( const ComplexType* x_val, ComplexType** deri );
      /*
       * Evaluates and differentiates the polynomials in the system
       * at the values in the array x_val.
       * On return is the function value of the system at x_val
       * and in deri is the Jacobian matrix evaluated at x_val.
       * The function assumes sufficient space is allocated for deri,
       * but allocates a new array of function values for the returned
       * values of the polynomial system at x_val.
       * The argument deri is the evaluated Jacobian matrix,
       * deri[idx] contains all evaluated partial derivatives of
       * the polynomial equation with index idx.
       * This function should not be used with multithreading.
       */

      void eval ( const ComplexType* x_val, ComplexType* f_val,
                  ComplexType** deri_val );
      /*
       * Evaluates and differentiates the polynomial system at x_val.
       * The function values are returned in f_val and all evaluated
       * partial derivatives of the polynomial equation with index idx
       * are returned in deri[idx].
       * If eval_base, then the powers of the common factor for each
       * equation are computed with memory allocation.
       * For general polynomial equations, this function should not be
       * applied in multithreaded runs.
       */

      void eval ( const ComplexType* x_val, ComplexType* f_val,
                  ComplexType** deri_val, ComplexType** deg_table );
      /*
       * Evaluates and differentiates the polynomial system at x_val.
       * The function values are returned in f_val and all evaluated
       * partial derivatives of the polynomial equation with index idx
       * are returned in deri[idx].
       * If eval_base, the powers of the common factor for each
       * equation are provided in deg_table.
       */

      void print();
      /*
       * Prints the polynomials in the system, using the print()
       * method defined in the class PolyEq.
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
      /*
       * Allocates space for the data attribute max_deg_base
       * and computes for each variable its highest power in
       * the system, with the method update_max_deg() in PolyEq.
       */

      ComplexType** eval_deg ( const ComplexType* x_val );
      /*
       * Computes the powers of the values for the coordinates in x_val,
       * using the exponents in max_deg_base.
       */

      void balance_eq ( const ComplexType* x_val );
      /*
       * Evaluates the polynomial system at the point with values
       * for its coordinates in x_val and then subtracts the result
       * of the evaluation from the constant of each polynomial.
       * After this operation, eval(x_val) should turn out zero.
       */
};

#include "polysys.tpp"

#endif
