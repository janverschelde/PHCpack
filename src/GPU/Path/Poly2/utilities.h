// utilities.h contains prototypes of templated utility functions
// the corresponding definitions are in the utilities.tpp file

#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

inline double read_number ( const char* number )
{
   return atof(number);
}

/* still todo: convert strings to double doubles and quad doubles ...
 
inline dd_real read_number ( const char* number )
{
   return dd_real(number);
}

inline qd_real read_number ( const char* number )
{
   return qd_real(number);
}
*/

void string_stop
 ( const string& mon_string, int& loc, const char* symbols,
   const int n_symbol, int l );
/*
 * Stop string at certain single symbol,
 * used by get_coef_complex. */

void read_until_line ( ifstream& myfile, string prefix );

void read_until_line ( ifstream& myfile, string* prefix, int n_prefix );

string get_number_string ( const string& mon_string, int& loc, int l );
/*
 * Read number string from string,
 * used by get_coef_complex.  */

template <class ComplexType, class RealType>
ComplexType get_complex_number ( ifstream& myfile );

template <class ComplexType>
ComplexType get_coef_complex ( const string& mon_string, int& loc );
/*
 * Read complex number from string.  */

template <class ComplexType>
ComplexType get_coef_real ( const string& mon_string, int& loc );

template <class RealType>
void print_matrix ( RealType** m, int row, int col );

template <class RealType>
void print_vector ( RealType* v, int n );

template <class RealType>
void print_result ( RealType* v, RealType** m, int dim, int n_eq );

void print_size ( size_t tmp_size_B );

double* rand_val ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex_frac ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex_one ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex_n ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex_unit ( int dim );

template <class ComplexType>
ComplexType* rand_val_complex_unit_n ( int dim );

string* var_list ( int dim, string v );

template <class RealType>
void print_number ( RealType v );

template <class ComplexType, class RealType>
void print_number_complex ( ComplexType v );

template <class ComplexType, class RealType>
void print_coef_complex ( ComplexType coef );

template <class ComplexType>
void cpu_speel0
 ( const ComplexType* x_val, unsigned short* pos, ComplexType* deri,
   const ComplexType& coef ); 

template <class ComplexType>
void cpu_speel
 ( const ComplexType* x_val, unsigned short* pos, ComplexType* deri,
   const ComplexType& coef );

template <class ComplexType>
void cpu_speel_with_base0
 ( const ComplexType* x_val, unsigned short* pos, unsigned short* exp,
   ComplexType* deri, const ComplexType& coef );

template <class ComplexType>
void cpu_speel_with_base
 ( const ComplexType* x_val, unsigned short* pos, unsigned short* exp,
   ComplexType* deri, const ComplexType& coef );

double time_interval1 ( struct timeval start, struct timeval end );

int log2ceil ( int n );

#include "utilities.tpp"

#endif
