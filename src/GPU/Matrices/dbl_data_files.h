/* The file dbl_data_files.h specifies functions to write matrices to file
 * and to read matrices from file, in double precision. */

#ifndef __dbl_data_files_h__
#define __dbl_data_files_h__

using namespace std;

void dbl_write_matrix ( string name, int dim, double **A );
/*
 * Writes the real matrix A of dimension dim to the file with name. */ 

void cmplx_write_matrix ( string name, int dim, double **Are, double **Aim );
/*
 * Writes the complex matrix A of dimension dim to the file with name.
 * All real parts are in Are, all imaginary parts in Aim.
 * The complex matrix on file is written as a sequence of
 * real and imaginary parts of each complex numbers in the matrix A. */

void dbl_read_matrix ( string name, int dim, double **A );
/*
 * Reads dim*dim real numbers of the file with name into A. */ 

void cmplx_read_matrix ( string name, int dim, double **Are, double **Aim );
/*
 * Reads dim*dim complex numbers of the file with name
 * into Are (the real parts) and Aim (the imaginary parts). */ 

#endif
