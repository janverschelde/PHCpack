// prototypes of code for execution on the host

#ifndef _gram_h
#define _gram_h

#include "DefineType.h"
#include "complex.h"
#include "complexH.h"

void print_gram_matrices ( complexH<T1> **g, complex<T> *g_h, int dim );
/*
 * DESCRIPTION :
 *   Prints the content of the two Gram matrices, g is computed by the host,
 *   and g_h is computed by the device, for visual inspection.
 *
 * ON ENTRY :
 *   g         matrix of dimension dim-by-dim of complex numbers;
 *   g_h       dim-by-dim matrix stored as one long array of size dim*dim;
 *   dim       dimension of the Gram matrices. */

void print_difference ( complexH<T1> **g, complex<T> *g_h, int dim );
/*
 * DESCRIPTION :
 *   Prints only the difference between the two Gram matrices,
 *   g is computed by the host and g_h is computed by the device,
 *   computed as the square root the sum of the differences between
 *   the components of the two matrices.
 *   Visual inspection is no longer feasible for larger dimensions.
 *
 * ON ENTRY :
 *   g         matrix of dimension dim-by-dim of complex numbers;
 *   g_h       dim-by-dim matrix stored as one long array of size dim*dim;
 *   dim       dimension of the Gram matrices. */

void CPU_inner_product
 ( complexH<T1>** v, int dim, int i, int j, complexH<T1>* ip );
/*
 * DESCRIPTION :
 *   Computes the complex conjugated inner product of v[i] with v[j].
 *
 * ON ENTRY :
 *   v         dim vectors of complex numbers of size dim;
 *   dim       dimension of the vectors in v;
 *   i         index to the first vector in v;
 *   j         index to the second vector in v.
 *
 * ON RETURN :
 *   ip        complex conjugated inner product of v[i] with v[j]. */

void CPU_gram ( complexH<T1>** v, complexH<T1>** g, int dim );
/*
 * DESCRIPTION :
 *   Computes the Gram matrix of a square matrix v.
 *
 * ON ENTRY :
 *   v         dim vectors of complex numbers of size dim;
 *   g         allocated space for a dim-by-dim matrix;
 *   dim       dimension of the vectors in v.
 *
 * ON RETURN :
 *   g         g[i][j] is the complex conjugated inner product
 *             of v[i] with v[j]. */

#endif
