/* file dc_roots.h contains the prototype for a function that returns
   approximations to all roots of a polynomial */

#include "dcmplx.h"

void roots ( int n, dcmplx p[n], double eps, int maxit, dcmplx r[n-1] );

/* DESCRIPTION : 
      The method of Durand-Kerner to approximate roots of a polynomial
      in one variable with complex coefficients.

   ON ENTRY :
      n        length of the coefficient vector, degree of p is n-1;
      p        p is the coefficient vector of a univariate polynomial,
               p[i] is the complex coefficient for the monomial x^i;
      eps      accuracy requirement on the roots;
      maxit    maximal number of allowed iterations.

   ON RETURN :
      r        approximations to the n-1 complex roots of p. */

void multiplicities ( int n, double tol, dcmplx r[n], int m[n] );
  
/* DESCRIPTION :
      A root has multiplicity m if there are m-1 roots close to it,
      in a cluster of radius the given tolerance tol.

   ON ENTRY :
      n        number of roots in r;
      tol      tolerance to decide whether two roots are equal;
      r        approximations for the roots.

   ON RETURN :
      m        multiplicities: m[i] is the multiplicity of r[i].  */

void multiple_roots ( int n, dcmplx p[n], double eps, int maxit, 
                      dcmplx r[n-1], double tol, int m[n-1] );

/* DESCRIPTION : 
      The method of Durand-Kerner to approximate roots of a polynomial
      in one variable with complex coefficients, followed by a Newton
      refinement on the (m-1)-th derivative, with m = multiplicity.

   ON ENTRY :
      n        length of the coefficient vector, degree of p is n-1;
      p        p is the coefficient vector of a univariate polynomial,
               p[i] is the complex coefficient for the monomial x^i;
      eps      accuracy requirement on the roots;
      maxit    maximal number of allowed iterations;
      tol      tolerance to decide whether two solutions are close
               to a multiple root.

   ON RETURN :
      r        approximations to the n-1 complex roots of p;
      m        m[i] is the multiplicity of the root in r[i]. */

int Newton ( int n, dcmplx p[n], dcmplx dp[n-1], dcmplx *z,
             double eps, int maxit );

/* DESCRIPTION :
      Newton's method to refine one root of a univariate polynomial.

   ON ENTRY :
      n        length of the coefficient vector, degree of p is n-1;
      p        p is the coefficient vector of a univariate polynomial,
               p[i] is the complex coefficient for the monomial x^i;
      dp       coefficient vector for the 1st derivative of p;
      z        pointer to an approximation for a root of p;
      eps      accuracy requirement: the iteration stops when the
               correction term becomes equal to or smaller than eps;
      maxit    maximal allowed number of iterations.

   ON RETURN :
      z        refinement of the given approximation for the root;
      Newton returns an integer with the following meaning:
        -1     failure to reach the accuracy requirement within 
               the maximal allowed number of iterations;
        otherwise, the integer on return is the number of iterations
                   that was needed to reach the accuracy requirement. */

dcmplx horner ( int n, dcmplx *p, dcmplx x );

/* DESCRIPTION :
      Horner's method to evaluate a polynomial at a point.

   ON ENTRY :
      n        length of the coefficient vector, degree of p is n-1;
      p        p[i] is the complex coefficient for the monomial x^i;
      x        complex number where p has to be evaluated.

   ON RETURN :
      p(x), the function value of p at the point x. */

dcmplx difference_product ( int n, int i, dcmplx z[n] );

/* DESCRIPTION :
      Auxiliary operation in the method of Durand-Kerner.

   ON ENTRY :
      n        number of elements in the vector z;
      i        an index between 0 and n-1;
      z        current approximations for the roots.

   ON RETURN :
      The product of all differences z[i]-z[j], j /= i. */

void derivative ( int n, dcmplx p[n], dcmplx dp[n-1] );

/* DESCRIPTION :
      Computes the coefficients of the derivative of the polynomial p.

   ON ENTRY :
      n        number of coefficients in the vector p;
      p        p[i] is the coefficient for the monomial x^i.

   ON RETURN :
      dp       coefficients of the derivative of p. */

void derivatives ( int n, int m, dcmplx p[n], dcmplx *dp[m] );

/* DESCRIPTION :
      Returns an array of the first m derivatives of the polynomial p.

   ON ENTRY :
      n        number of coefficients in the vector p;
      m        number of derivatives wanted;
      p        p[i] is the coefficient for the monomial x^i.

   ON RETURN :
      dp       an array with the first m derivatives of p, dp[i] is 
               the coefficient vector of the (i+1)-th derivative. */

void write_derivatives ( int n, int m, dcmplx *dp[m] );

/* DESCRIPTION :
      Writes the first m derivatives of a polynomial to screen.

   ON ENTRY :
      n        length of the coefficient vector of the original polynomial;
      m        number of derivatives;
      dp       coefficient vector of the derivatives. */
  
