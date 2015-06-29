#ifndef __FAMILY_H__
#define __FAMILY_H__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "DefineType_Host.h"

using namespace std;

/** @file 
@brief Generator of families of Polynomial system functions*/

/// Generate random doubles for n variables
/**
@param dim dimension or size of array
@return random double array of size dim
*/
double* rand_val(int dim);

CT* rand_val_complex(int dim);

CT* rand_val_complex_one(int dim);

CT* rand_val_complex_n(int dim);

CT* rand_val_complex_frac(int dim);

CT* rand_val_complex_unit(int dim);

CT* rand_val_complex_unit_n(int dim);

/// Generate strings of variables
/**
Variable string = x + dim
@param dim dimension or size of array
@return stirng array of size dim
*/
string* x_var(string x, int dim);

/// Write a cyclic-n funtion into file.
/**
Generate a cyclic polynomial system of dimension dim,
and write into the file, f_name.
@return void
@param f_name string of file name
@param dim integer of dimension
*/
void file_cyclic(string f_name, int dim);

/// Generate a cyclic-n function as a string array.
/**
Generate a cyclic polynomial system of dimension dim into a string array.
@param f_name string of file name
@param dim integer of dimension
@return string array of equations in cyclic-n
*/
string* string_cyclic(int dim);

string* string_two(int dim);


#endif
