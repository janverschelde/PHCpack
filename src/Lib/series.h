/* The file series.h contains prototypes to the function to compute power
 * series solutions of polynomial systems by Newton's method.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SERIES_H__
#define __SERIES_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int standard_Newton_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision.  The solution series are stored in 
 *   the standard systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   nbr     number of Newton steps to be done on each solution;
 *   maxdeg  the maximal degree of the series;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int dobldobl_Newton_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision.  The solution series are stored in the
 *   double double systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   maxdeg  the maximal degree of the series;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int quaddobl_Newton_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision.  The solution series are stored in the
 *   quad double systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   maxdeg  the maximal degree of the series;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int standard_Newton_power_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the standard systems pool
 *   the leading terms of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision.  The solution series are stored in the
 *   standard systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   maxdeg  the maximal degree of the series;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int dobldobl_Newton_power_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in double double precision, and in the dobldobl systems pool
 *   the leading terms of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision.  The solution series are stored in the
 *   double double systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   maxdeg  the maximal degree of the series;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int quaddobl_Newton_power_series ( int idx, int maxdeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in quad double precision, and in the quad double systems pool
 *   the leading terms of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision.  The solution series are stored in the 
 *   quad double systems pool.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   maxdeg  the maximal degree of the series;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int standard_Pade_approximant
 ( int idx, int numdeg, int dendeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision.  The solution series are stored in 
 *   the standard systems pool.  The series lead to Pade approximants.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   numdeg  degree of the numerator;
 *   dendeb  degree of the denominator;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int dobldobl_Pade_approximant
 ( int idx, int numdeg, int dendeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in double double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision.  The solution series are stored in 
 *   the dobldobl systems pool.  The series lead to Pade approximants.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   numdeg  degree of the numerator;
 *   dendeb  degree of the denominator;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int quaddobl_Pade_approximant
 ( int idx, int numdeg, int dendeg, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in quad double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision.  The solution series are stored in 
 *   the quaddobl systems pool.  The series lead to Pade approximants.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   numdeg  degree of the numerator;
 *   dendeb  degree of the denominator;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

#endif
