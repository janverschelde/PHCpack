/* file schubert.h contains prototypes to the numerical Schubert calculus
   available in PHCpack */

int Pieri_root_count ( int m, int p, int q, int *r );
/*
 * Returns in r the number of p-plane producing curves
 * of degree q the meet m-planes at prescribed points. */

int localization_poset ( int m, int p, int q, int *nc, char *ps );
/*
 * Returns in ps, a string of nc characters, the string representation
 * of the localization poset for the Pieri root count for (m,p,q). */

int Pieri_polynomial_system
 ( int m, int p, int q, int nc, char *A, int is_real );
/*
 * Makes the polynomial system defined by the planes in A.
 *
 * ON ENTRY :
 *   m        dimension of the input planes;
 *   p        dimension of the output planes;
 *   q        degree of the solution maps;
 *   nc       number of characters in the string A;
 *   A        m*p + q*(m+p) random complex input m-planes,
 *            the real and imaginary parts are separated by a space;
 *   is_real  if == 1, then the coefficients of A are real,
 *            if == 0, then the coefficients of A are complex.
 * 
 * The system container contains the polynomial system 
 * defined by the input planes in A. */

int run_Pieri_homotopies
 ( int m, int p, int q, int nc, int *r, char *A, char *pts );
/*
 * Runs the Pieri homotopies for (m,p,q) dimensions on random input data.
 *
 * ON ENTRY :
 *   m        dimension of the input planes;
 *   p        dimension of the output planes;
 *   q        degree of the solution maps;
 *   nc       number of characters in the string A;
 *   A        m*p + q*(m+p) random complex input m-planes,
 *            the real and imaginary parts are separated by a space;
 *   pts      m*p + q*(m+p) random complex interpolation points,
 *            only if q > 0.
 *
 * ON RETURN :
 *   r        the combinatorial Pieri root count.
 * 
 * The system container contains the polynomial system solved
 * and the solutions are in the solution container. */

int real_osculating_planes
 ( int m, int p, int q, int *nc, char *s, char *planes );
/*
 * Returns in planes (a string of nc characters), the string representation
 * of n real m-planes in d-space osculating a rational normal curve
 * at the n points in s, where n = m*p + q*(m+p) and d = m+p. */
