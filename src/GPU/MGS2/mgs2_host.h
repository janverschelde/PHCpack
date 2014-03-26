// Prototypes of code for execution on the host fall into three categories:
// 1. print and copy functions;
// 2. check normality, orthogonality, decomposition, and solution
//    for data computed by the host, as well as by the device;
// 3. routines for orthonormalization, a QR decomposition, and
//    QR followed by back substitution to solve a linear system.

#ifndef _mgs2_h
#define _mgs2_h

#include "DefineType.h"
#include "complex.h"
#include "complexH.h"

void print_complex_vector ( complex<T> *v, int dim );
/*
 * DESCRIPTION :
 *   Prints the vector of dimension dim. */

void print_real_vector ( T *v, int dim );
/*
 * DESCRIPTION :
 *   Prints the vector of dimension dim. */

void print_complex_vector ( complexH<T1> *v, int dim );
/*
 * DESCRIPTION :
 *   Prints the vector of dimension dim. */

void print_real_vector2 ( T1 *v, int dim );
/*
 * DESCRIPTION :
 *   Prints the vector of dimension dim. */

void print_complex_matrix ( complexH<T1> **v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the content of the rows-by-cols matrix v. */

void print_real_matrix ( T1 **v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the content of the rows-by-cols matrix v. */

void print_complex_matrices
 ( complexH<T1>** v, complex<T>* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the content of the two matrices, v is computed by the host,
 *   and v_h is computed by the device, for visual inspection.
 *
 * ON ENTRY :
 *   v        matrix of dimension rows-by-cols of complex numbers,
 *            the matrix is stored column wise, that is v[i] contains
 *            the entries of the i-th column of v;
 *   v_h      rows-by-cols matrix stored as one long array of size rows*cols;
 *   rows     number of rows of the matrices v and v_h;
 *   cols     number of colums of the matrix v and v_h. */

void print_real_matrices ( T1** v, T* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the content of the two matrices, v is computed by the host,
 *   and v_h is computed by the device, for visual inspection.
 *
 * ON ENTRY :
 *   v        matrix of dimension rows-by-cols of real numbers,
 *            the matrix is stored column wise, that is v[i] contains
 *            the entries of the i-th column of v;
 *   v_h      rows-by-cols matrix stored as one long array of size rows*cols;
 *   rows     number of rows of the matrices v and v_h;
 *   cols     number of colums of the matrix v and v_h. */

void complex_print_difference
 ( complexH<T1> **v, complex<T>* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints only the difference between the two matrices,
 *   v is computed by the host and v_h is computed by the device,
 *   computed as the square root the sum of the differences between
 *   the components of the two matrices.
 *   Visual inspection is no longer feasible for larger dimensions.
 *
 * ON ENTRY :
 *   v        matrix of dimension rows-by-cols of complex numbers;
 *            the matrix is stored column wise, that is v[i] contains
 *            the entries of the i-th column of v;
 *   v_h      rows-by-cols matrix stored as one long array of size rows*cols;
 *   rows     number of rows of the matrices v and v_h;
 *   cols     number of colums of the matrix v and v_h. */

void real_print_difference
 ( T1 **v, T *v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints only the difference between the two matrices,
 *   v is computed by the host and v_h is computed by the device,
 *   computed as the square root the sum of the differences between
 *   the components of the two matrices.
 *   Visual inspection is no longer feasible for larger dimensions.
 *
 * ON ENTRY :
 *   v        matrix of dimension rows-by-cols of real numbers;
 *            the matrix is stored column wise, that is v[i] contains
 *            the entries of the i-th column of v;
 *   v_h      rows-by-cols matrix stored as one long array of size rows*cols;
 *   rows     number of rows of the matrices v and v_h;
 *   cols     number of colums of the matrix v and v_h. */

void copy_complex_matrices
 ( complexH<T1>** vfrom, complexH<T1>** vto, int rows, int cols );
/*
 * DESCRIPTION :
 *   Copies the matrix vfrom to the matrix vto.
 *
 * REQUIRED :
 *   Both matrices vfrom and vto are row-by-cols matrices.
 *
 * ON ENTRY :
 *   vfrom    row-by-cols matrix of complex numbers;
 *   vto      space allocated for a complex row-by-cols matrix.
 *
 * ON RETURN :
 *   vto      row-by-cols matrix with same content as vfrom. */

void copy_real_matrices ( T1** vfrom, T1** vto, int rows, int cols );
/*
 * DESCRIPTION :
 *   Copies the matrix vfrom to the matrix vto.
 *
 * REQUIRED :
 *   Both matrices vfrom and vto are row-by-cols matrices.
 *
 * ON ENTRY :
 *   vfrom    row-by-cols matrix of real numbers;
 *   vto      space allocated for a real row-by-cols matrix.
 *
 * ON RETURN :
 *   vto      row-by-cols matrix with same content as vfrom. */

void complex_checkGPUnormal ( complex<T>* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v_h have normal one,
 *   printing the accumulated differences between one and the computed norms.
 *
 * ON ENTRY :
 *   v_h      rows*cols complex numbers of matrix stored column wise;
 *   rows     number of rows of each column in v_h;
 *   cols     number of columns in v_h. */

void real_checkGPUnormal ( T* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v_h have normal one,
 *   printing the accumulated differences between one and the computed norms.
 *
 * ON ENTRY :
 *   v_h      rows*cols real numbers of matrix stored column wise;
 *   rows     number of rows of each column in v_h;
 *   cols     number of columns in v_h. */

void complex_checkCPUnormal ( complexH<T1>** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v have normal one,
 *   printing the accumulated differences between one and the computed norms.
 *
 * ON ENTRY :
 *   v        rows*cols complex numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void real_checkCPUnormal ( T1** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v have normal one,
 *   printing the accumulated differences between one and the computed norms.
 *
 * ON ENTRY :
 *   v        rows*cols real numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void complex_checkGPUorthogonal ( complex<T>* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v_h are mutually orthogonal,
 *   printing the sum of the computed inner products between the columns.
 *
 * ON ENTRY :
 *   v        rows*cols complex numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void real_checkGPUorthogonal ( T* v_h, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v_h are mutually orthogonal,
 *   printing the sum of the computed inner products between the columns.
 *
 * ON ENTRY :
 *   v        rows*cols real numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void complex_checkCPUorthogonal ( complexH<T1>** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v are mutually orthogonal,
 *   printing the sum of the computed inner products between the columns.
 *
 * ON ENTRY :
 *   v        rows*cols complex numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void real_checkCPUorthogonal ( T1** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the columns in the matrix v are mutually orthogonal,
 *   printing the sum of the computed inner products between the columns.
 *
 * ON ENTRY :
 *   v        rows*cols real numbers of matrix stored column wise;
 *   rows     number of rows of each column in v;
 *   cols     number of columns in v. */

void complex_checkGPUdecomposition
 ( complexH<T1>** A, complex<T>* Q, complex<T>* R,
   int dimR, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the matrix A equals Q^T*R, prints the sum of the
 *   differences between all components of A and Q^T*R.
 *
 * ON ENTRY :
 *   A        rows-by-cols matrix of complex numbers stored column wise;
 *   Q        orthonormal basis for the column span of A;
 *   R        multipliers used in the orthonormalization of A;
 *   dimR     number of multipliers equals cols*(cols+1)/2;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void real_checkGPUdecomposition
 ( T1** A, T* Q, T* R, int dimR, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the matrix A equals Q^T*R, prints the sum of the
 *   differences between all components of A and Q^T*R.
 *
 * ON ENTRY :
 *   A        rows-by-cols matrix of real numbers stored column wise;
 *   Q        orthonormal basis for the column span of A;
 *   R        multipliers used in the orthonormalization of A;
 *   dimR     number of multipliers equals cols*(cols+1)/2;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void complex_checkCPUdecomposition
 ( complexH<T1>** A, complexH<T1>** Q, complexH<T1>** R, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the matrix A equals Q^T*R, prints the sum of the
 *   differences between all components of A and Q^T*R.
 *
 * ON ENTRY :
 *   A        rows*cols complex number of matrix stored column wise;
 *   Q        orthonormal basis for the column span of A;
 *   R        multipliers used in the orthonormalization of A;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void real_checkCPUdecomposition
 ( T1** A, T1** Q, T1** R, int rows, int cols );
/*
 * DESCRIPTION :
 *   Checks whether the matrix A equals Q^T*R, prints the sum of the
 *   differences between all components of A and Q^T*R.
 *
 * ON ENTRY :
 *   A        rows*cols real number of matrix stored column wise;
 *   Q        orthonormal basis for the column span of A;
 *   R        multipliers used in the orthonormalization of A;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void complex_checkGPUsolution
 ( complexH<T1>** A, complex<T>* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the residual of the solution x of the linear system
 *   defined by the columns in A, where the last column of A is
 *   the right hand side vector of the linear system.
 *
 * ON ENTRY :
 *   A        rows*cols complex number of matrix stored column wise;
 *   x        vector of cols-1 complex numbers with a solution;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void real_checkGPUsolution ( T1** A, T* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the residual of the solution x of the linear system
 *   defined by the columns in A, where the last column of A is
 *   the right hand side vector of the linear system.
 *
 * ON ENTRY :
 *   A        rows*cols real number of matrix stored column wise;
 *   x        vector of cols-1 real numbers with a solution;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void complex_checkCPUsolution
 ( complexH<T1>** A, complexH<T1>* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the residual of the solution x of the linear system
 *   defined by the columns in A, where the last column of A is
 *   the right hand side vector of the linear system.
 *
 * ON ENTRY :
 *   A        rows*cols complex numbers of matrix stored column wise;
 *   x        vector of cols-1 complex numbers with a solution;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void real_checkCPUsolution ( T1** A, T1* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Prints the residual of the solution x of the linear system
 *   defined by the columns in A, where the last column of A is
 *   the right hand side vector of the linear system.
 *
 * ON ENTRY :
 *   A        rows*cols real numbers of matrix stored column wise;
 *   x        vector of cols-1 real numbers with a solution;
 *   rows     number of rows of each column in A;
 *   cols     number of columns in A. */

void CPU_complex_normalize_and_reduce
 ( complexH<T1>** v, int rows, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Normalizes the pivot column of v and reduces all later columns in v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise,
 *             with an orthonormal basis in columns 0 to pivot-1;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v;
 *   pivot     current column in v to be normalized
 *             and used for reduction of the later columns.
 *
 * ON RETURN :
 *   v         orthonormal basis in columns 0 to pivot. */

void CPU_real_normalize_and_reduce
 ( T1** v, int rows, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Normalizes the pivot column of v and reduces all later columns in v.
 *
 * ON ENTRY :
 *   v         matrix of real numbers stored column wise,
 *             with an orthonormal basis in columns 0 to pivot-1;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v;
 *   pivot     current column in v to be normalized
 *             and used for reduction of the later columns.
 *
 * ON RETURN :
 *   v         orthonormal basis in columns 0 to pivot. */

void CPU_complex_QR_normalize_and_reduce
 ( complexH<T1>** v, complexH<T1>** R, int rows, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Normalizes the pivot column of v and reduces all later columns in v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise,
 *             with an orthonormal basis in columns 0 to pivot-1;
 *   R         matrix allocated for a cols-by-cols matrix;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v;
 *   pivot     current column in v to be normalized
 *             and used for reduction of the later columns.
 *
 * ON RETURN :
 *   v         orthonormal basis in columns 0 to pivot;
 *   R         the pivot column of R contains the multipliers:
 *             R[pivot][pivot] is the 2-norm of the pivot column,
 *             R[pivot][k] for k > pivot contains the inner product
 *             of the pivot column with the k-th column of v,
 *             so R is stored as a lower triangular matrix. */

void CPU_real_QR_normalize_and_reduce
 ( T1** v, T1** R, int rows, int cols, int pivot );
/*
 * DESCRIPTION :
 *   Normalizes the pivot column of v and reduces all later columns in v.
 *
 * ON ENTRY :
 *   v         matrix of real numbers stored column wise,
 *             with an orthonormal basis in columns 0 to pivot-1;
 *   R         matrix allocated for a cols-by-cols matrix;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v;
 *   pivot     current column in v to be normalized
 *             and used for reduction of the later columns.
 *
 * ON RETURN :
 *   v         orthonormal basis in columns 0 to pivot;
 *   R         the pivot column of R contains the multipliers:
 *             R[pivot][pivot] is the 2-norm of the pivot column,
 *             R[pivot][k] for k > pivot contains the inner product
 *             of the pivot column with the k-th column of v,
 *             so R is stored as a lower triangular matrix. */

void CPU_complex_backsubstitution
 ( complexH<T1>** U, complexH<T1>* rhs, complexH<T1>* x, int dim );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = rhs, of dimension dim.
 *
 * ON ENTRY :
 *   U         square upper triangular matrix of dimension dim;
 *   rhs       right-hand size vector;
 *   x         memory allocated for dim complex numbers;
 *   dim       dimension of the linear system.
 *
 * ON RETURN :
 *   x         solution of the system U*x = rhs. */

void CPU_real_backsubstitution ( T1** U, T1* rhs, T1* x, int dim );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = rhs, of dimension dim.
 *
 * ON ENTRY :
 *   U         square upper triangular matrix of dimension dim;
 *   rhs       right-hand size vector;
 *   x         memory allocated for dim real numbers;
 *   dim       dimension of the linear system.
 *
 * ON RETURN :
 *   x         solution of the system U*x = rhs. */

void CPU_complex_mgs2 ( complexH<T1>** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the complex matrix v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal. */

void CPU_real_mgs2 ( T1** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the real matrix v.
 *
 * ON ENTRY :
 *   v         matrix of real numbers stored column wise;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal. */

void CPU_complex_mgs2qr
 ( complexH<T1>** v, complexH<T1>** R, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the modified Gram-Schmidt method 
 *   to compute the QR decomposition of the matrix v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise;
 *   R         matrix allocated for a rows-by-cols matrix;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal;
 *   R         contains the multipliers, Q*A = R,
 *             where A is the v on input and Q the v on output. */

void CPU_real_mgs2qr ( T1** v, T1** R, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the modified Gram-Schmidt method 
 *   to compute the QR decomposition of the matrix v.
 *
 * ON ENTRY :
 *   v         matrix of real numbers stored column wise;
 *   R         matrix allocated for a rows-by-cols matrix;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal;
 *   R         contains the multipliers, Q*A = R,
 *             where A is the v on input and Q the v on output. */

void CPU_complex_mgs2qrls 
 ( complexH<T1>** v, complexH<T1>** R, complexH<T1>* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the matrix v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise;
 *   R         matrix allocated for a rows-by-cols matrix;
 *   x         space for cols-1 complex numbers;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal;
 *   R         contains the multipliers, Q*A = R,
 *             where A is the v on input and Q the v on outputi;
 *   x         solution to the system defined by the first cols-1 columns
 *             in v with the last column of v as right hand size vector. */

void CPU_real_mgs2qrls ( T1** v, T1** R, T1* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the matrix v.
 *
 * ON ENTRY :
 *   v         matrix of real numbers stored column wise;
 *   R         matrix allocated for a rows-by-cols matrix;
 *   x         space for cols-1 real numbers;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal;
 *   R         contains the multipliers, Q*A = R,
 *             where A is the v on input and Q the v on outputi;
 *   x         solution to the system defined by the first cols-1 columns
 *             in v with the last column of v as right hand size vector. */

#endif
