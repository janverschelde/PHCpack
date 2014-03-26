// Prototypes to evaluate and differentiate the polynomial equations
// arising from a discretization of the Chandrasekhar H-equation,
// following the formulation in a paper by Laureano Gonzalez-Vega: 
// "Some examples on problem solving by using the
//  symbolic viewpoint when dealing with polynomial systems of equations".
// In: "Computer Algebra in Science and Engineering", pages 102-116,
// edited by J. Fleischer, J. Grabmeier, F.W. Hehl and W. Kuechlin.
// World Scientific Publishing, 1995.

#include <iostream>
#include "DefineType.h"
#include "complexH.h"

using namespace std;

void complex_chandra_evaluate
 ( complexH<T1> c, int dim, complexH<T1>* x, complexH<T1>* y );
/*
 * DESCRIPTION :
 *   Returns the value of the chandra system for dimension dim,
 *   for the parameter c, evaluated at x with results in y.
 *
 * REQUIRED : space has been allocated in y for dim complex numbers. */

void real_chandra_evaluate ( T1 c, int dim, T1* x, T1* y );
/*
 * DESCRIPTION :
 *   Returns the value of the chandra system for dimension dim,
 *   for the parameter c, evaluated at x with results in y.
 *
 * REQUIRED : space has been allocated in y for dim real numbers. */

void complex_chandra_differentiate
 ( complexH<T1> c, int dim, int i, int j, complexH<T1>* x, complexH<T1>* y );
/*
 * DESCRIPTION :
 *   Returns the value of the derivative with respect to the j-th variable
 *   of the i-th equation in the chandra system for dimension dim for the
 *   parameter c, evaluated at x with result in y.
 *   This computes the (i,j)-th component of the Jacobian matrix.
 *
 * REQUIRED : space has been allocated in y for a complex number. */

void real_chandra_differentiate ( T1 c, int dim, int i, int j, T1* x, T1* y );
/*
 * DESCRIPTION :
 *   Returns the value of the derivative with respect to the j-th variable
 *   of the i-th equation in the chandra system for dimension dim for the
 *   parameter c, evaluated at x with result in y.
 *   This computes the (i,j)-th component of the Jacobian matrix.
 *
 * REQUIRED : space has been allocated in y for a real number. */

void complex_chandra_evaluate_and_differentiate
 ( complexH<T1> c, int dim, complexH<T1>* x, complexH<T1>** v );
/*
 * DESCRIPTION :
 *   Returns in v the linear system for running Newton's method.
 *
 * ON ENTRY :
 *   c         parameter of the chandra system;
 *   dim       dimension of the system;
 *   x         dim complex numbers;
 *   v         space allocated for dim+1 columns for dimension dim;
 *             in a column oriented matrix: v[0] contains the first column
 *             v[1] the second column, etc.
 *
 * ON RETURN :
 *   v         v[i] contains the i-th column of the Jacobian matrix,
 *             evaluated at x, for i ranging from 0 to dim-1
 *             v[dim] contains the value of the chandra system at x,
 *             with flipped sign. */

void real_chandra_evaluate_and_differentiate ( T1 c, int dim, T1* x, T1** v );
/*
 * DESCRIPTION :
 *   Returns in v the linear system for running Newton's method.
 *
 * ON ENTRY :
 *   c         parameter of the chandra system;
 *   dim       dimension of the system;
 *   x         dim real numbers;
 *   v         space allocated for dim+1 columns for dimension dim;
 *             in a column oriented matrix: v[0] contains the first column
 *             v[1] the second column, etc.
 *
 * ON RETURN :
 *   v         v[i] contains the i-th column of the Jacobian matrix,
 *             evaluated at x, for i ranging from 0 to dim-1
 *             v[dim] contains the value of the chandra system at x,
 *             with flipped sign. */
