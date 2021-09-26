/* The file dbl_data_files.h specifies functions to write matrices to file
 * and to read matrices from file, in double precision. */

#ifndef __dbl_data_files_h__
#define __dbl_data_files_h__

using namespace std;

void dbl_write_matrix ( string name, int dim, double **A );
/*
 * Writes the matrix A of dimension dim to the file with name. */ 

void dbl_reads_matrix ( string name, int dim, double **A );
/*
 * Reads dim*dim numbers of the file with name into A. */ 

#endif
