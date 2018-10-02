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

// Monomial class
/*
   Read the monomial from the string.
   PolyMon is used to construct polynomial equation, PolyEq,
   which is used to construct polynomial system, PolySys.
   @sa PolyEq PolySys
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
      // Makes a monomial for n variables, with exponents in d,
      // and real and imaginary parts of the coefficient in c.
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

      ~PolyMon()
      {
         if(pos != NULL) delete[] pos;
         if(exp != NULL) delete[] exp;
      }

      void read 
         ( const string& eq_string, VarDict& pos_dict,
           int start, int end, ComplexType coef );
      // Reads a monomial from a certain part in the equation string.
      /*
         Coefficient is given.
         @param eq_string equation string
         @param pos_dict the dictionary of variables and their positions
         @param start start index of the monomial in equation string
         @param end end index of the monomial in equation string
         @param coef0 the coefficient of the monomail
       */

      void read ( const string& mon_string, VarDict& pos_dict,
                  ComplexType coef );
      // Reads a monomial from a monomial string.
      /*
         Coefficient is given.
         @param mon_string monomial string
         @param pos_dict the dictionary of variables and their positions
         @param coef0 the coefficient of the monomail
       */

      void read ( const string& mon_string, VarDict& pos_dict );
      // Reads a monomial from a monomial string.
      /*
         Read coefficient first and then variables.
         @param mon_string monomial string
         @param pos_dict the dictionary of variables and their positions
       */

      ComplexType speel ( const ComplexType* x_val, ComplexType* deri );
      /*
         Applies algorithm of Speelpenning to evaluate the monomial and all
         its derivatives.  The first n_var positions are deri are filled.
         This function applies only when n_base == 0, that is:
         when no variables appear with power 2 or higher.
         The function assumes the monomial contains at least one variable
         with an exponent equal to one.
       */
      /*
         First compute forward product, 
         then compute backward and cross product together.
         @param[in] x_val values of variables
         @param[out] deri array of derivatives
         @return value of product of variables
       */

      ComplexType speel_with_base
       ( const ComplexType* x_val, ComplexType* deri, ComplexType base );

      ComplexType eval_base
       ( const ComplexType* x_val, ComplexType** deg_table );
      // Compute base of monomial
      /*
         @param[in] x_val values of variables
         @return value of base
       */

      ComplexType eval ( const ComplexType* x_val );
      // The monomial is evaluated with the straightforward algorithm.
      /*
         @param x_val array of variables' values
         @return ComplexType monomial value
       */

      ComplexType eval ( const ComplexType* x_val, ComplexType* deri );
      // Evaluate monomial and its derivatives
      /*
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return monomial value
         @sa eval_base() speel()
       */

      ComplexType eval
       ( const ComplexType* x_val, ComplexType* deri,
         ComplexType** deg_table );

      void print ( const string* pos_var );
      // prints a monomial
      /*
         Variable symbols are read from pos_dict.
         @param pos_dict the dictionary of variables and their positions
       */

      int memory_size ( int factor_size );

      void print_level(); // prints the level structure

      int job_number_block ( int start_level );

      void update_max_deg ( int* max_deg );

      void update_base();
      // Auxiliary operation for when reading a monomial.
};

#include "polymon.tpp"

#endif
