// Prototypes of code for execution on the host fall into three categories:
// 1. print and copy functions;
// 2. check normality, orthogonality, decomposition, and solution
//    for data computed by the host, as well as by the device;
// 3. routines for orthonormalization, a QR decomposition, and
//    QR followed by back substitution to solve a linear system.
// The matrices are stored in a column wise fashion, that is:
// the column index is the leading index.

#ifndef __MGS2_H__
#define __MGS2_H__

template <class ComplexType, class RealType>
void CPU_normalize_and_reduce
 ( ComplexType** v, int rows, int cols, int pivot );
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

template <class ComplexType, class RealType>
void CPU_QR_normalize_and_reduce
 ( ComplexType** v, ComplexType** R, int rows, int cols, int pivot );
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

template <class ComplexType>
void CPU_backsubstitution
 ( ComplexType** U, ComplexType* rhs, ComplexType* x, int dim );
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

template <class ComplexType, class RealType>
void CPU_mgs2 ( ComplexType** v, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the matrix v.
 *
 * ON ENTRY :
 *   v         matrix of complex numbers stored column wise;
 *   rows      number of rows in the matrix v,
 *             equals the size of v[i], for i in range 0..cols-1;
 *   cols      number of columns in the matrix v.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal. */

template <class ComplexType, class RealType>
void CPU_mgs2qr ( ComplexType** v, ComplexType** R, int rows, int cols );
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

template <class ComplexType, class RealType>
void CPU_mgs2qrls
 ( ComplexType** v, ComplexType** R, ComplexType* x, int rows, int cols,
   ComplexType* rhs );
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
 *   cols      number of columns in the matrix v;
 *   rhs       workspace allocated for the right hand side vector,
 *             with as many elements as the value of cols.
 *
 * ON RETURN :
 *   v         the vectors in the columns of v are orthonormal;
 *   R         contains the multipliers, Q*A = R,
 *             where A is the v on input and Q the v on outputi;
 *   x         solution to the system defined by the first cols-1 columns
 *             in v with the last column of v as right hand size vector;
 *   rhs       updated workspace vector. */

template <class ComplexType, class RealType>
void CPU_mgs2qrls
 ( ComplexType** v, ComplexType** R, ComplexType* x, int rows, int cols );
/*
 * DESCRIPTION :
 *   Performs the Gram-Schmidt method on the matrix v.
 *   Allocates and deallocates space for the right hand side vector rhs.
 *   Wraps the previous CPU_mgs2qrls.
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

#include "mgs_host.tpp"

#endif
