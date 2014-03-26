// Prototypes for the kernels for the modified Gram-Schmidt method
// on a sequence of complex or real vectors are of three kinds:
// 1. Functions that have 'small' in their name are for dimensions where the
//    problem data fits into shared memory.  For double, double double,
//    and quad double precision, the maximum small dimension is
//    respectively 256, 128, and 64.  For runs on small data, the
//    block size (number of threads per block) equals the dimension.
// 2. Functions that have 'large' in their name are for dimensions that can
//    be up to 32 times larger than 256, 128, and 64 for double,
//    double double, and quad double precision.  For runs on large data,
//    the block size remains limited to 256, 128, or 64, for the double
//    double double, and quad double precision respectively.
// 3. The functions GPU_mgs2, GPU_mgs2qr, and GPU_mgs2qrls are responsible
//    for the data transfer and launching of the kernel functions,
//    respectively for an orthonormalization, a QR decompsition, and
//    the solving of a linear system following a QR decomposition.

#ifndef __MGS2_KERNELS__
#define __MGS2_KERNELS__

#include <iostream>
#include <cmath>
#include <gqd_type.h>
#include "DefineType.h"
#include "complex.h"
#include "vector_types.h"
#include "vector_functions.h"

__global__ void complex_small_normalize_and_reduce
 ( complex<T>* v, int rows, int rowsLog2, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v
 *   and to reduce the later columns, for small dimensions.
 *
 * REQUIRED :
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   v         rows*cols complex numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the pivot column in v that will be normalized,
 *             the block index points to the vector that will be reduced.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized and every column
 *             indexed by the block is reduced with this vector. */

__global__ void real_small_normalize_and_reduce
 ( T* v, int rows, int rowsLog2, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v
 *   and to reduce the later columns, for small dimensions.
 *
 * REQUIRED :
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   v         rows*cols real numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the pivot column in v that will be normalized,
 *             the block index points to the vector that will be reduced.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized and every column
 *             indexed by the block is reduced with this vector. */

__global__ void complex_large_normalize
 ( complex<T>* v, int rows, int pivot, int BS, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v,
 *   for large dimension, done by as many blocks as the value of rnd.
 *
 * ON ENTRY :
 *   v         rows*cols complex numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   pivot     index to the current pivot column in v,
 *             this is the column in v that will be normalized;
 *   BS        the block size or the number of threads in the block;
 *   pivnorm   norm of the column pivot.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized. */

__global__ void real_large_normalize
 ( T* v, int rows, int pivot, int BS, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v,
 *   for large dimension, done by as many blocks as the value of rnd.
 *
 * ON ENTRY :
 *   v         rows*cols real numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   pivot     index to the current pivot column in v,
 *             this is the column in v that will be normalized;
 *   BS        the block size or the number of threads in the block;
 *   pivnorm   norm of the column pivot.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized. */

__global__ void complex_large_normalize_and_reduce
 ( complex<T>* v, int rows, int rowsLog2, int cols, int pivot,
   int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the previous pivot column of v
 *   and to reduce the later columns, for larger dimensions.
 *
 * ON ENTRY :
 *   v         rows*cols complex numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the current pivot column in v,
 *             the block index points to the vector that will be reduced;
 *   rnd       the number of rounds or the multiplier for the number of rows;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction;
 *   pivnorm   norm of the column pivot-1, in case pivot > 0.
 *
 * ON RETURN :
 *   v         matrix with previous pivot column normalized and columns
 *             indexed by the block index are reduced with this vector;
 *   pivnorm   norm of the current pivot column. */

__global__ void real_large_normalize_and_reduce
 ( T* v, int rows, int rowsLog2, int cols, int pivot,
   int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the previous pivot column of v
 *   and to reduce the later columns, for larger dimensions.
 *
 * ON ENTRY :
 *   v         rows*cols real numbers for cols vectors of dimension rows;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the current pivot column in v,
 *             the block index points to the vector that will be reduced;
 *   rnd       the number of rounds or the multiplier for the number of rows;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction;
 *   pivnorm   norm of the column pivot-1, in case pivot > 0.
 *
 * ON RETURN :
 *   v         matrix with previous pivot column normalized and columns
 *             indexed by the block index are reduced with this vector;
 *   pivnorm   norm of the current pivot column. */

__global__ void complex_small_QR_normalize_and_reduce
 ( complex<T>* v, complex<T>* R,
   int dimR, int rows, int rowsLog2, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v, to reduce
 *   the later columns, and to store multipliers, for small dimensions.
 *
 * REQUIRED :
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   v         rows*cols complex numbers for cols vectors of dimension rows;
 *   R         space allocated for dimR complex numbers;
 *   dimR      equals cols*(cols+1)/2, the dimension of R;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the pivot column in v that will be normalized,
 *             the block index points to the vector that will be reduced.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized and every column
 *             indexed by the block is reduced with this vector;
 *   R         multipliers stored for coalesced access in back substitution. */

__global__ void real_small_QR_normalize_and_reduce
 ( T* v, T* R, int dimR, int rows, int rowsLog2, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the pivot column of v, to reduce
 *   the later columns, and to store multipliers, for small dimensions.
 *
 * REQUIRED :
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   v         rows*cols real numbers for cols vectors of dimension rows;
 *   R         space allocated for dimR real numbers;
 *   dimR      equals cols*(cols+1)/2, the dimension of R;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the pivot column in v that will be normalized,
 *             the block index points to the vector that will be reduced.
 *
 * ON RETURN :
 *   v         matrix with pivot column normalized and every column
 *             indexed by the block is reduced with this vector;
 *   R         multipliers stored for coalesced access in back substitution. */

__global__ void complex_large_QR_normalize_and_reduce
 ( complex<T>* v, complex<T>* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the previous pivot column of v
 *   to reduce the later columns, and storing the multipliers for
 *   a QR decomposition of v, for larger dimensions.
 *
 * ON ENTRY :
 *   v         rows*cols complex numbers for cols vectors of dimension rows;
 *   R         space allocated for dimR complex numbers;
 *   dimR      equals cols*(cols+1)/2, the dimension of R;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the current pivot column in v,
 *             the block index points to the vector that will be reduced;
 *   rnd       the number of rounds or the multiplier for the number of rows;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction;
 *   pivnorm   norm of the column pivot-1, in case pivot > 0.
 *
 * ON RETURN :
 *   v         matrix with previous pivot column normalized and columns
 *             indexed by the block index are reduced with this vector;
 *   R         multipliers stored for coalesced access in back substitution;
 *   pivnorm   norm of the current pivot column. */

__global__ void real_large_QR_normalize_and_reduce
 ( T* v, T* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm );
/*
 * DESCRIPTION :
 *   Kernel function to normalize the previous pivot column of v
 *   to reduce the later columns, and storing the multipliers for
 *   a QR decomposition of v, for larger dimensions.
 *
 * ON ENTRY :
 *   v         rows*cols real numbers for cols vectors of dimension rows;
 *   R         space allocated for dimR real numbers;
 *   dimR      equals cols*(cols+1)/2, the dimension of R;
 *   rows      number of rows of every vector in v;
 *   rowsLog2  equals ceil(log2((double) rows), used in sum reduction.
 *   cols      number of columns in the matrix v;
 *   pivot     index to the current pivot column in v,
 *             the block index points to the vector that will be reduced;
 *   rnd       the number of rounds or the multiplier for the number of rows;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction;
 *   pivnorm   norm of the column pivot-1, in case pivot > 0.
 *
 * ON RETURN :
 *   v         matrix with previous pivot column normalized and columns
 *             indexed by the block index are reduced with this vector;
 *   R         multipliers stored for coalesced access in back substitution;
 *   pivnorm   norm of the current pivot column. */

__global__ void complex_small_backsubstitution
 ( complex<T>* R, complex<T>* x, int dim );
/*
 * DESCRIPTION :
 *   Performs a back substitution on the linear system defined
 *   by the data in R, for small dimensions.
 *
 * REQUIRED : 
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 * 
 * ON ENTRY :
 *    R        output of the QR decomposition with as last column
 *             the reduced right hand side vector of the linear system;
 *    x        space allocated for dim complex numbers;
 *    dim      dimension of the system.
 *
 * ON RETURN :
 *    x        the solution to the linear system defined by R. */

__global__ void real_small_backsubstitution ( T* R, T* x, int dim );
/*
 * DESCRIPTION :
 *   Performs a back substitution on the linear system defined
 *   by the data in R, for small dimensions.
 *
 * REQUIRED : 
 *   The number of rows equals the block size,
 *   which is the number of threads in the block.
 * 
 * ON ENTRY :
 *    R        output of the QR decomposition with as last column
 *             the reduced right hand side vector of the linear system;
 *    x        space allocated for dim real numbers;
 *    dim      dimension of the system.
 *
 * ON RETURN :
 *    x        the solution to the linear system defined by R. */

__global__ void complex_large_backsubstitution
 ( complex<T>* R, complex<T>* x, int dim, int rnd, int pivot, int BS  );
/*
 * DESCRIPTION :
 *   Performs a back substitution on the linear system defined by R,
 *   for any dimension, running in as many rounds as the value of rnd.
 *
 * ON ENTRY :
 *    R        output of the QR decomposition with as last column
 *             the reduced right hand side vector of the linear system;
 *    x        space allocated for dim complex numbers;
 *    dim      dimension of the system;
 *    rnd      number of rounds;
 *    pivot    the pivot indicates the triangular subsystem of size BS
 *             that will be solved by block 0, and the components of the
 *             solution are indexed by pivot + BS + j, for j = 0, .., BS-1;
 *    BS       the number of threads in one block.
 *
 * ON RETURN :
 *    x        the solution to the linear system defined by R. */

__global__ void real_large_backsubstitution
 ( T* R, T* x, int dim, int rnd, int pivot, int BS  );
/*
 * DESCRIPTION :
 *   Performs a back substitution on the linear system defined by R,
 *   for any dimension, running in as many rounds as the value of rnd.
 *
 * ON ENTRY :
 *    R        output of the QR decomposition with as last column
 *             the reduced right hand side vector of the linear system;
 *    x        space allocated for dim real numbers;
 *    dim      dimension of the system;
 *    rnd      number of rounds;
 *    pivot    the pivot indicates the triangular subsystem of size BS
 *             that will be solved by block 0, and the components of the
 *             solution are indexed by pivot + BS + j, for j = 0, .., BS-1;
 *    BS       the number of threads in one block.
 *
 * ON RETURN :
 *    x        the solution to the linear system defined by R. */

__global__ void complex_chandra_evaluate_and_differentiate
 ( complex<T>* v, complex<T>* R, complex<T>* x, int dimR, int rows,
   int rowsLog2, int cols, int rnd, int rndLog2, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates the discretized Chandrasekhar equation
 *   at x into the v and R matrices for use in Newton's method. */

void GPU_complex_mgs2
 ( complex<T>* v_h, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols complex numbers for cols vectors of size rows;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h. */

void GPU_real_mgs2 ( T* v_h, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols real numbers for cols vectors of size rows;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h. */

void GPU_complex_mgs2qr
 ( complex<T>* v_h, complex<T>* R_h,
   int dimR, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method
 *   to compute the QR decomposition of the vectors stored in v_h,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols complex numbers for cols vectors of size rows;
 *   R_h       space allocated for the dimR multipliers;
 *   dimR      equals cols*(cols+1)/2;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h,
 *   R_h       matrix of multipliers so A = Q*R. */

void GPU_real_mgs2qr
 ( T* v_h, T* R_h, int dimR, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method
 *   to compute the QR decomposition of the vectors stored in v_h,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols real numbers for cols vectors of size rows;
 *   R_h       space allocated for the dimR multipliers;
 *   dimR      equals cols*(cols+1)/2;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h,
 *   R_h       matrix of multipliers so A = Q*R. */

void GPU_complex_mgs2qrls
 ( complex<T>* v_h, complex<T>* R_h, complex<T> *x_h,
   int dimR, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method
 *   to compute the QR decomposition of the vectors stored in v_h,
 *   applies back substitution to solve the linear system,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols complex numbers for cols vectors of size rows,
 *             the last column is the right hand size vector of the linear
 *             system with coefficient matrix in the first cols-1 columns;
 *   R_h       space allocated for the dimR multipliers;
 *   x_h       space allocated as many complex numbers as rows;
 *   dimR      equals cols*(cols+1)/2;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h,
 *   R_h       matrix of multipliers so A = Q*R.
 *   x_h       the least-squares solution to the linear system defined
 *             by the columns of v_h as given on input. */

void GPU_real_mgs2qrls
 ( T* v_h, T* R_h, T *x_h, int dimR, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, runs the modified Gram-Schmidt method
 *   to compute the QR decomposition of the vectors stored in v_h,
 *   applies back substitution to solve the linear system,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       rows*cols real numbers for cols vectors of size rows,
 *             the last column is the right hand size vector of the linear
 *             system with coefficient matrix in the first cols-1 columns;
 *   R_h       space allocated for the dimR multipliers;
 *   x_h       space allocated as many real numbers as rows;
 *   dimR      equals cols*(cols+1)/2;
 *   rows      number of rows in every vector of v;
 *   cols      number of columns in v;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h       orthornormal basis for the column space of v_h,
 *   R_h       matrix of multipliers so A = Q*R.
 *   x_h       the least-squares solution to the linear system defined
 *             by the columns of v_h as given on input. */

void GPU_complex_Newton_chandra
 ( complex<T>* v_h, complex<T>* R_h, complex<T> *x_h,
   int dimR, int rows, int cols, int freq, int BS );
/*
 * DESCRIPTION :
 *   This function performs as many Newton steps as the value of freq
 *   on the discretization of the H-Chadrasekhar equation.
 *   All operations in Newton's method are performed by the GPU. */

#endif
