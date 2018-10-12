// The file test_utils.h contains the prototypes of test utility functions
// to test the operations on monomials as defined by polymon, and
// to test the operations on polynomials as defined by polyeq.
// The utility functions are applied in test_polymon and test_polyeq.

#ifndef __TEST_UTILS_H__
#define __TEST_UTILS_H__

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "complexH.h"
#include "polymon.h"
#include "polyeq.h"

template <class ComplexType, class RealType>
ComplexType random_complex ();
/*
 * Returns a random complex number. */

template <class ComplexType, class RealType>
void random_point ( int dim, ComplexType* coordinates );
/*
 * Generates the coordinates of a random complex point.
 * Space must be allocated for the dim many coordinates. */

void random_exponents ( int dim, int expmax, int* exponents );
/*
 * Generates dim random positive integers as exponents in [0, expmax],
 * uniformly distributed.
 * Space must be allocated in exponents for dim integers. */

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> random_monomial ( int dim, int expmax );
/*
 * Returns a random monomial, with a random complex coefficient,
 * for dim variables, with exponents in [0, expmax]. */

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> prompted_monomial ( int dim );
/*
 * Prompts the user for a coefficient and dim exponents
 * to define a monomial which is returned. */

template <class ComplexType, class RealType>
ComplexType plain_eval ( PolyMon<ComplexType,RealType>& m, ComplexType *x );
/*
 * Applies the straightforward algorithm to evaluate m at x. */

template <class ComplexType, class RealType>
void plain_diff
 ( PolyMon<ComplexType,RealType>& m, ComplexType *x, ComplexType *deri );
/*
 * Applies the straightforward algorithm to differentiate m at x.
 * The values of the derivatives are returned in deri in their order
 * of appearance in the monomial, at positions 0 to m.n_var. */

template <class ComplexType, class RealType>
void print_data ( PolyMon<ComplexType,RealType>& monomial );
/*
 * Prints the content of the data stored in monomial. */

template <class ComplexType, class RealType>
PolyEq<ComplexType,RealType> random_polynomial
 ( int dim, int nbterms, int expmax );
/*
 * Returns a random polynomial in a space with dim variables
 * and with as many terms as the value of nbterms.
 * The largest exponent is defined by the value of expmax. */

template <class ComplexType, class RealType>
void write_monomial ( int dim, PolyMon<ComplexType,RealType>& m );
/*
 * Writes the information stored in the monomial m in tableau format, 
 * first the real and imaginary part of the coefficient,
 * followed by all dim exponents, where dim is the ambient dimension.
 * There is no newline written at the end.
 */

template <class ComplexType, class RealType>
void write_polynomial ( PolyEq<ComplexType,RealType>& p );
/*
 * Writes the terms in p in tableau style format. */

template <class ComplexType, class RealType>
ComplexType plain_eval
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x );
/*
 * Applies the straightforward algorithm to evaluate p at x. */

template <class ComplexType, class RealType>
void plain_diff
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x, ComplexType *deri );
/*
 * Applies the straightforward algorithm to differentiate p at x.
 * The function returns in deri the values of all derivatives. */

template <class ComplexType, class RealType>
void powertable
 ( int dim, int* maxdeg, ComplexType* point, ComplexType** powers );
/*
 * Computes the table with the powers of the variables,
 * at the dim coordinates of the given point,
 * given the maximal degrees for each variable from 0 to dim-1.
 * The function assumes sufficient space for powers has been allocated.
 * Row idx in powers[idx] contains the powers of the variable idx,
 * starting at the value point[idx]. */

#include "test_utils.tpp"

#endif
