// The file cyclic_columns.h specifies function to make the column
// representation of the cyclic n-roots system.

#ifndef __cyclic_columns_h__
#define __cyclic_columns_h__

void make_cyclic_variables ( int nbrcol, int dim, int **nvr );
/*
 * DESCRIPTION :
 *   Defines the number of variables in each monomial
 *   in the column representation of the cyclic n-roots system.
 *
 * ON ENTRY :
 *   nbrcol  is the number of columns, where nbrcol <= dim,
 *           if nbrcol == dim, then the cyclic n-roots is complete;
 *   dim     number of equations and total number of variables;
 *   nvr     space for dim columns and dim variables in each column.
 *
 * ON RETURN :
 *   nvr     nvr[i][j] is the number of variables of the j-th monomial
 *           in the i-th column. */

void make_cyclic_columns ( int nbrcol, int dim, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Defines the column representation of the cyclic n-roots system.
 *
 * ON ENTRY :
 *   nbrcol  is the number of columns, complete if nbrcol == dim;
 *   dim     number of equations and total number of variables;
 *   nvr     nvr[[i][j] is the number of variables of the j-th monomial
 *           in the i-th column;
 *   idx     space for the variables that appear in every monomial
 *           in every column.
 *
 * ON RETURN :
 *   idx     idx[i][j][k] is the index of the k-th variable which appears
 *           in the j-th monomial of the i-th column. */

void write_cyclic_columns ( int nbrcol, int dim, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Writes the exponents of the column representation of
 *   the cyclic n-roots systems,
 *   for a number of columns equals to nbrcol. */

#endif
