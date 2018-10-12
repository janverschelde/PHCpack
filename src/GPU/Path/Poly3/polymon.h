// The file polymon.h contains prototypes of templated data types for
// monomials with complex coefficients of different precisions.
// The corresponding definitions of the functions are in polymon.tpp.

#ifndef __POLYMON_H__
#define __POLYMON_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "dict.h"
#include "utilities.h"

using namespace std;

/*
 * The class PolyMon defines the representation of a monomial,
 * provides methods to read a monomial from string, and
 * methods to efficiently evaluate and differentiate a monomial
 * at a given point.
 * The coefficient type is complex and the templated construction
 * allows for different precisions.
 */

template <class ComplexType, class RealType>
class PolyMon
{
   public:

      ComplexType coef;  // coefficient of the monomial
      int n_var;         // number of variables with exponent > 0
      int dim;           // ambient dimension, total number of variables
      int* pos;          // array of size n_var with positions of variables
      int* exp;          // array of size n_var with exponents of variables

      int n_base;        // number of variables with exponent > 1
      int* pos_base;     // array of size n_base with positions of variables
      int* exp_base;     // array of size n_base with exponents of variables
      int* exp_tbl_base; // stores values of exponents minus 2

      PolyMon() // makes the zero constant
      {
         coef = ComplexType(0.0,0.0);
         n_var = 0;
         dim = 0;
         pos = NULL;
         exp = NULL;
         n_base = 0;
         pos_base = NULL;
         exp_base = NULL;
         exp_tbl_base = NULL;
      }

      PolyMon ( int dim ) // the zero constant in dimension dim
      {
         coef = ComplexType(0.0,0.0);
         n_var = 0;
         this -> dim = dim;
         pos = NULL;
         exp = NULL;
         n_base = 0;
         pos_base = NULL;
         exp_base = NULL;
         exp_tbl_base = NULL;
      }

      PolyMon ( int n, int* d, RealType* c )
      /*
       * Makes a monomial for n variables, with exponents in d,
       * and real and imaginary parts of the coefficient in c.
       */
      {
         coef = ComplexType(c[0],c[1]);
         dim = n;
         n_var = 0;
         n_base = 0;
         for(int var_idx=0; var_idx<dim; var_idx++)
         {
            if(d[var_idx]!=0)
            {
               n_var++;
               if(d[var_idx]>1) n_base++;
            }
         }
         pos = new int[n_var];
         exp = new int[n_var];
         pos_base = new int[n_base];
         exp_base = new int[n_base];
         exp_tbl_base = new int[n_base];
         int pos_idx = 0;
         int base_idx = 0;
         for(int var_idx=0; var_idx<dim; var_idx++)
         {
            if(d[var_idx]!=0)
            {
               pos[pos_idx] = var_idx;
               exp[pos_idx] = d[var_idx];
               pos_idx++;
               if(d[var_idx]>1)
               {
                  pos_base[base_idx] = var_idx;
                  exp_base[base_idx] = d[var_idx];
                  exp_tbl_base[base_idx] = d[var_idx]-2;
                  base_idx++;
               }
            }
         }
      }

      ~PolyMon() // The destructor frees the allocated space.
      {
         if(pos != NULL) delete[] pos;
         if(exp != NULL) delete[] exp;
         if(pos_base != NULL) delete[] pos_base;
         if(exp_base != NULL) delete[] exp_base;
         if(exp_base != NULL) delete[] exp_tbl_base;
      }

      void read ( const string& eq_string, VarDict& pos_dict,
                  int start, int end, ComplexType coef, int verbose=0 );
      /*
       * Reads the string from position start to end for the exponents
       * of a monomial written in symbolic representation in the string.
       * The dictionary pos_dict defines the positions of the variables.
       * The coefficient of the monomial is given in coef.
       * If verbose, then extra output is written to screen.
       */ 

      void read ( const string& mon_string, VarDict& pos_dict,
                  ComplexType coef, int verbose=0 );
      /*
       * Given the coefficient coef of a monomial, reads the string
       * which contains the monomial in symbolic format for the exponents.
       * The dictionary pos_dict define the position of the variables.
       * If verbose, then extra output is written to screen.
       */

      void read ( const string& mon_string, VarDict& pos_dict, int verbose=0 );
      /*
       * Reads a monomial from a string which contains a monomial
       * in symbolic representation.  The dictionary pos_dict defines
       * the positions of the variables.
       * If verbose, then extra output is written to screen.
       */

      ComplexType speel ( const ComplexType* x_val, ComplexType* deri );
      /*
       * Applies algorithm of Speelpenning to evaluate the monomial and all
       * its derivatives.  The first n_var positions are deri are filled.
       * This function applies only when n_base == 0, that is:
       * when no variables appear with power 2 or higher.
       * The function assumes the monomial contains at least one variable
       * with an exponent equal to one.
       */

      ComplexType speel_with_base
       ( const ComplexType* x_val, ComplexType* deri, ComplexType base );
      /*
       * Applies the Speelpenning algorithm to evaluate and differentiate
       * the monomial at the point defined by x_val.
       * The value of the monomial at x_val is returned.
       * The derivatives are returned in deri.
       * The value base is the common factor defined by all higher degree
       * powers of the variables in the monomial.
       */

      ComplexType eval_base
       ( const ComplexType* x_val, ComplexType** deg_table );
      /*
       * Multiplies the coefficient of the monomial with the powers
       * of the variables in the double array deg_table.
       * While the values in x_val are not used, the first entries
       * in deg_table should equal to the coordinates in xval.
       * Returns the value of the coefficient of the monomial
       * times the value of the common factor of the monomial.
       */

      ComplexType eval ( const ComplexType* x_val );
      /*
       * The monomial is evaluated with the straightforward algorithm
       * at x_val, which must be an array of range 0..dim-1.
       */

      ComplexType eval ( const ComplexType* x_val, ComplexType* deri );
      /*
       * Wraps the method speel() to evaluate and differentiate the
       * monomial at the point x_val, with derivatives returned in deri.
       */

      ComplexType eval
       ( const ComplexType* x_val, ComplexType* deri,
         ComplexType** deg_table );
      /*
       * Wraps the method speel_with_base() to evaluate and differentiate
       * the monomial at the point x_val, with powers of the variables
       * given in deg_table.  The derivatives are returned in deri.
       */

      void print ( const string* pos_var );
      /*
       * Prints the monomial, given the positions of the variables
       * defined in the string pos_var.
       */

      void print_tableau ( int dim );
      /*
       * Writes the information stored in the monomial m in tableau format, 
       * first the real and imaginary part of the coefficient,
       * followed by all dim exponents, where dim is the ambient dimension.
       * There is no newline written at the end.
       */

      int memory_size ( int factor_size );

      void print_level(); // prints the level structure

      int job_number_block ( int start_level );

      void update_max_deg ( int* max_deg );
      /*
       * Updates the array of integers in max_deg with the value
       * of the largest exponent, relative to its position.
       */

      void update_base();
      // Auxiliary operation for when reading a monomial.
};

#include "polymon.tpp"

#endif
