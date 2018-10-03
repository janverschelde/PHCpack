// The file test_utils.h contains the prototypes of test utility functions
// to test the operations on monomials as defined by polymon.

#ifndef __TEST_UTILS_H__
#define __TEST_UTILS_H__

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "complexH.h"
#include "polymon.h"

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

#include "test_utils.tpp"

#endif
