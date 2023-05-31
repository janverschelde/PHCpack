// The file cyclic_indices.h specifies function to define the indices 
// of the monomials in the cyclic n-roots system.

#ifndef __cyclic_indices_h__
#define __cyclic_indices_h__

void make_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Makes the indices of the cyclic n-roots system.
 *
 * ON ENTRY :
 *   dim      the number of polynomials in the system.
 *   nbr      space for dim numbers;
 *   nvr      space for dim pointers;
 *   idx      space for dim pointers.
 *
 * ON RETURN :
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial.    */

void write_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Writes the indices of a polynomial system.
 *
 * ON ENTRY :
 *   dim      number of polynomials in the system;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial.    */

#endif
