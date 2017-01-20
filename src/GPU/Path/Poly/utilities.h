#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>

#include "DefineType_Host.h"

using namespace std;

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

CT get_complex_number ( ifstream& myfile );

CT get_coef_complex ( const string& mon_string, int& loc );
/*
 * Read complex number from string.  */

CT get_coef_real ( const string& mon_string, int& loc );

template <class T1>
void print_matrix ( T1** m, int row, int col );

template <class T1>
void print_vector ( T1* v, int n );

template <class T1>
void print_result ( T1* v, T1** m, int dim, int n_eq );

void print_size ( size_t tmp_size_B );

double* rand_val ( int dim );

CT* rand_val_complex ( int dim );

CT* rand_val_complex_frac ( int dim );

CT* rand_val_complex_one ( int dim );

CT* rand_val_complex_n ( int dim );

CT* rand_val_complex_unit ( int dim );

CT* rand_val_complex_unit_n ( int dim );

string* var_list ( int dim, string v );

void print_number ( T1 v );

void print_number_complex ( CT v );

void print_coef_complex ( CT coef );

void cpu_speel0
 ( const CT* x_val, unsigned short* pos, CT* deri, const CT& coef );

void cpu_speel
 ( const CT* x_val, unsigned short* pos, CT* deri, const CT& coef );

void cpu_speel_with_base0
 ( const CT* x_val, unsigned short* pos, unsigned short* exp, CT* deri,
   const CT& coef );

void cpu_speel_with_base
 ( const CT* x_val, unsigned short* pos, unsigned short* exp, CT* deri,
   const CT& coef );

double time_interval1 ( struct timeval start, struct timeval end );

int log2ceil ( int n );

#endif
