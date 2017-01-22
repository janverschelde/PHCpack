#ifndef __POLY_H__
#define __POLY_H__

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

using namespace std;

// Monomial class
/*
Read the monomial from the string.
PolyMon is used to construct polynomial equaion, PolyEq,
which is used to construct polynomial system, PolySys.
@sa PolyEq PolySys
*/

class PolyMon
{
   public:

      CT coef; // coefficient of the monomial
      int n_var; // derivative array of the monomial
      int dim; // position array of variables in the monomial
      int* pos;
      int* exp; // exponent array of variables in the monomial

      int n_base;
      int* pos_base;
      int* exp_base;
      int* exp_tbl_base;

      PolyMon()
      {
         coef = CT(0.0,0.0);
         n_var = 0;
         dim = 0;
         pos = NULL;
         exp = NULL;
         n_base = 0;
         pos_base = NULL;
         exp_base = NULL;
         exp_tbl_base = NULL;
      }

      PolyMon ( int dim )
      {
         coef = CT(0.0,0.0);
         n_var = 0;
         this -> dim = dim;
         pos = NULL;
         exp = NULL;
         n_base = 0;
         pos_base = NULL;
         exp_base = NULL;
         exp_tbl_base = NULL;
      }

      PolyMon ( int n, int* d, T1* c )
      {
         coef = CT(c[0],c[1]);
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
         delete[] pos;
         delete[] exp;
      }

      void read 
         ( const string& eq_string, VarDict& pos_dict,
           int start, int end, CT coef );
      // Reads a monomial from a certain part in the equation string.
      /*
         Coefficient is given.
         @param eq_string equation string
         @param pos_dict the dictionary of variables and their positions
         @param start start index of the monomial in equation string
         @param end end index of the monomial in equation string
         @param coef0 the coefficient of the monomail
       */

      void read ( const string& mon_string, VarDict& pos_dict, CT coef );
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

      CT speel ( const CT* x_val, CT* deri );
      // Use speel expanding to compute derivatives
      /*
         First compute forward product, 
         then compute backward and cross product together.
         @param[in] x_val values of variables
         @param[out] deri array of derivatives
         @return value of product of variables
       */

      CT speel_with_base ( const CT* x_val, CT* deri, CT base );

      CT eval_base ( const CT* x_val, CT** deg_table );
      // Compute base of monomial
      /*
         @param[in] x_val values of variables
         @return value of base
       */

      CT eval ( const CT* x_val ); // Evaluate monomial
      /*
         @param x_val array of variables' values
         @return CT monomial value
       */

      CT eval ( const CT* x_val, CT* deri );
      // Evaluate monomial and its derivatives
      /*
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return monomial value
         @sa eval_base() speel()
       */

      CT eval ( const CT* x_val, CT* deri, CT** deg_table );

      void print ( const string* pos_var );
      // prints a monomial
      /*
         Variable symbols are read from pos_dict.
         @param pos_dict the dictionary of variables and their positions
       */

      int memory_size(int factor_size);

      void print_level();
      // print level structure

      int job_number_block ( int start_level );

      void update_max_deg(int* max_deg);

      void update_base();
};

/// Polynomial equation class
/*
Read polynomial equation from string,
and build array of monomials, PolyMon.
PolyEq is used to construct polynomial system, PolySys.
@sa PolyMon PolySys
*/
class PolyEq
{
   public:

      CT constant; // Constant of equation
      int n_mon;   // Number of monomials
      int dim;     // Number of monomials
      vector<PolyMon*> mon; // Monomial array of the equation

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
          for(vector<PolyMon*>::iterator it=mon.begin();
              it<mon.end(); it++)
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

      CT eval ( const CT* x_val );
      // Evaluate equation
      /*
         @param x_val CT array of variables' values
         @return CT equation value
       */
	
      CT eval ( const CT* x_val, CT* deri );
      // Evaluate equation and its derivatives
      /*
         Evaluate monomials and their derivative by PolyMon::eval(),
         and add them up.
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return equation value
       */

      CT eval(const CT* x_val, CT* deri, CT** deg_table);

      int memory_size ( int factor_size );

      void print_level();

      int workspace_size_block ( int start_level, int factor_size );

      int job_number_block ( int start_level );

      void update_max_deg ( int* max_deg );
};

// Polynomial system class
/*
Read polynomial system from string or string arrays,
and build array of polynomial equations, PolyEq.
@sa PolyMon, PolyEq.
*/
class PolySys
{
   public:

      int n_eq; // Number of Equations
      int dim;  // Dimension
      string* pos_var;
      vector<PolyEq*> eq; // Array of polynomial equations
	
      PolyEq* eq_space;

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
         delete max_deg_base;
         // cout << "sys destructed" << endl;
      }

      void read ( const string* sys_string, int n_eq, VarDict& pos_dict );
      // Read polynomial system from string array
      /*
        Read equation from string array, sys_string, using PolyEq::read().
        @param sys_string polynomial system string
        @param n number of equations or number of strings
        in string array sys_string
        @param pos_dict the dictionary of variables and their positions
        @sa PolyEq::read()
       */

      void read ( const string& sys_string, VarDict& pos_dict );
      // Read polynomial system from string
      /*
        First split the string by ';',
        then use PolyEq::read() to read each equation.
        @param sys_string polynomial system string
        @param pos_dict the dictionary of variables and their positions
        @sa PolyEq::read()
       */

      void read_file ( const string& file_name );
      // Read a polynomial system from file.
      /*
         Read dimension on the first line and then system.
         @param filename file name of polynomial system.
         File requirement: first line is dimension
         and the equations are separated by ';' symbols.
         @param pos_dict position dictionary for variables.
         If it is empty, the positions will be created by the order 
         of variables appearance first time in the file.
         @sa read_file(ifstream& myfile, Dict& pos_dict)
       */

      void read_file ( ifstream& myfile, VarDict& pos_dict );
      // Read polynomial system from file
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

      CT* eval ( const CT* x_val );
      // Evaluate polynomial system
      /*
         @param x_val CT array of variables' values
         @return CT array of system value
      */

      CT* eval ( const CT* x_val, CT** deri );
      // Evaluate equation and its derivative
      /*
         Evaluate monomials and their derivatives by PolyMon::eval(),
         and add them up.
         @param[in] x_val array of variables' values
         @param[out] deri array of derivatives
         @return equation value
       */

      void eval ( const CT* x_val, CT* f_val,  CT** deri_val );

      void print();
      // Print polynomial system
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

      CT** eval_deg ( const CT* x_val );

      void balance_eq ( const CT* x_val );
};

class PolySysHom
{
   public:

      PolySys* start_sys;
      PolySys* target_sys;
      int dim;

      PolySysHom ( PolySys* start_sys, PolySys* target_sys )
      {
         if(start_sys->dim != target_sys->dim)
         {
            std::cout << "start system and end system " << std::endl;
         }
         else
         {
            this->start_sys = start_sys;
            this->target_sys = target_sys;
            dim = start_sys->dim;
         }
      }

      void print()
      {
         std::cout << "Start System : " << std::endl;
         start_sys->print();
         std::cout << "Target System : " << std::endl;
         target_sys->print();
      }
};

class int_idx
{
   public:

      int eq_idx;
      int mon_idx;

      int_idx ( int i, int j )
      {
         eq_idx = i;
         mon_idx = j;
      }
};

#endif
