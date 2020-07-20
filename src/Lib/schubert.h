/* The file schubert.h contains prototypes to the numerical Schubert calculus
 * available in PHCpack.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SCHUBERT_H__
#define __SCHUBERT_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int Pieri_root_count ( int m, int p, int q, int *r );
/*
 * Returns in r the number of p-plane producing curves
 * of degree q the meet m-planes at prescribed points. */

int resolve_Schubert_conditions
 ( int n, int k, int c, int *brackets, int verbose, int *r );
/*
 * In n-space for k-planes subject to c general Schubert intersection
 * conditions as defined by the brackets, the root count r is computed.
 *
 * ON ENTRY :
 *   n        the ambient dimension, where the k-planes live;
 *   k        the dimension of the solution planes;
 *   c        the number of intersection conditions;
 *   brackets is an array of k*c integer numbers with the conditions
 *            imposed by the flags on the k-plane;
 *   verbose  when 0, no intermediate output is written,
 *            when 1, then the resolution is dispayed on screen.
 *
 * ON RETURN :
 *   r        the formal root count as the number of k-planes
 *            for conditions imposed by the brackets for general flags. */

int standard_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare, 
   int nbchar, char *filename, int *r, double *flags );
/*
 * In n-space for k-planes subject to c general Schubert intersection
 * conditions as defined by the brackets, the root count r is computed
 * and the Littlewood-Richardson homotopies compute as many k-planes
 * as the value of r for generic flags, in standard double precision.
 * Upon return, the polynomial system that was solved is in the systems
 * container and its solutions are in the solutions container.
 *
 * ON ENTRY :
 *   n        the ambient dimension, where the k-planes live;
 *   k        the dimension of the solution planes;
 *   c        the number of intersection conditions;
 *   brackets is an array of k*c integer numbers with the conditions
 *            imposed by the flags on the k-plane;
 *   verify   if 0, no diagnostic verification is done,
 *            if 1, diagnostic verification output is written to file;
 *   verbose  if 0, no intermediate output is written,
 *            if 1, then the resolution is dispayed on screen;
 *   minrep   if 0, all minors are used in the problem formulation,
 *            if 1, a more efficient problem formulation is used;
 *   tosquare if 0, then the overdetermined homotopies are not squared,
 *            if 1, then squaring is applied to the overdetermined systems;
 *   nbchar   number of characters in the string filename;
 *   filename is the name of the output file.
 *
 * ON RETURN :
 *   r        the formal root count as the number of k-planes
 *            for conditions imposed by the brackets for general flags;
 *   flags    coefficients of the general flags. */

int dobldobl_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare,
   int nbchar, char *filename, int *r, double *flags );
/*
 * In n-space for k-planes subject to c general Schubert intersection
 * conditions as defined by the brackets, the root count r is computed
 * and the Littlewood-Richardson homotopies compute as many k-planes
 * as the value of r for generic flags, in double double precision.
 * Upon return, the polynomial system that was solved is in the systems
 * container and its solutions are in the solutions container.
 *
 * ON ENTRY :
 *   n        the ambient dimension, where the k-planes live;
 *   k        the dimension of the solution planes;
 *   c        the number of intersection conditions;
 *   brackets is an array of k*c integer numbers with the conditions
 *            imposed by the flags on the k-plane;
 *   verify   if 0, no diagnostic verification is done,
 *            if 1, diagnostic verification output is written to file;
 *   verbose  if 0, no intermediate output is written,
 *            if 1, then the resolution is dispayed on screen;
 *   minrep   if 0, all minors are used in the problem formulation,
 *            if 1, a more efficient problem formulation is used;
 *   tosquare if 0, then the overdetermined homotopies are not squared,
 *            if 1, then squaring is applied to the overdetermined systems;
 *   nbchar   number of characters in the string filename;
 *   filename is the name of the output file.
 *
 * ON RETURN :
 *   r        the formal root count as the number of k-planes
 *            for conditions imposed by the brackets for general flags;
 *   flags    coefficients of the general flags. */

int quaddobl_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare,
   int nbchar, char *filename, int *r, double *flags );
/*
 * In n-space for k-planes subject to c general Schubert intersection
 * conditions as defined by the brackets, the root count r is computed
 * and the Littlewood-Richardson homotopies compute as many k-planes
 * as the value of r for generic flags, in quad double precision.
 * Upon return, the polynomial system that was solved is in the systems
 * container and its solutions are in the solutions container.
 *
 * ON ENTRY :
 *   n        the ambient dimension, where the k-planes live;
 *   k        the dimension of the solution planes;
 *   c        the number of intersection conditions;
 *   brackets is an array of k*c integer numbers with the conditions
 *            imposed by the flags on the k-plane;
 *   verify   if 0, no diagnostic verification is done,
 *            if 1, diagnostic verification output is written to file;
 *   verbose  if 0, no intermediate output is written,
 *            if 1, then the resolution is dispayed on screen;
 *   minrep   if 0, all minors are used in the problem formulation,
 *            if 1, a more efficient problem formulation is used;
 *   tosquare if 0, then the overdetermined homotopies are not squared,
 *            if 1, then squaring is applied to the overdetermined systems;
 *   nbchar   number of characters in the string filename;
 *   filename is the name of the output file.
 *
 * ON RETURN :
 *   r        the formal root count as the number of k-planes
 *            for conditions imposed by the brackets for general flags;
 *   flags    coefficients of the general flags. */

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

#endif
