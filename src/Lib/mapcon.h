/* This file "mapcon.h" contains the prototypes of the operations
 * on the monomial maps container in PHCpack.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __MAPCON_H__
#define __MAPCON_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int mapcon_solve_system ( int puretopdim );
/* 
 * DESCRIPTION :
 *   Solves the binomial system stored in the Laurent systems container.
 *   If puredim equals 1, then only the pure top dimensional solution sets
 *   will be computed.
 *   If puredim equals 0, all solution sets are returned.*/

int mapcon_write_maps ( void );
/* 
 * DESCRIPTION :
 *   Writes the maps stored in the container to screen. */

int mapcon_clear_maps ( void );
/* 
 * DESCRIPTION :
 *   Destroys the maps stored in the container. */

int mapcon_top_dimension ( int *dim );
/* 
 * DESCRIPTION :
 *   Returns in dim the top dimension of the maps in the container. */

int mapcon_number_of_maps ( int dim, int *nbmaps );
/* 
 * DESCRIPTION :
 *   Returns in nbmaps the number of maps of dimension dim. */

int mapcon_degree_of_map ( int dim, int ind, int *deg );
/* 
 * DESCRIPTION :
 *   Returns in deg the degree of the map of dimension dim and index ind. */

int mapcon_coefficients_of_map ( int dim, int ind, int nvr, double *cff );
/* 
 * DESCRIPTION :
 *   Returns in cff the coefficients of the map of dimension dim
 *   and index ind, where nvr is the number of variables.
 *
 * REQUIRED : cff is an array of doubles with 2*nvr allocated space
 *   for the real and imaginary parts of the coefficients. */

int mapcon_exponents_of_map ( int dim, int ind, int nvr, int *exp );
/* 
 * DESCRIPTION :
 *   Returns in exp the exponents of the map of dimension dim
 *   and index ind, where nvr is the number of variables.
 *
 * REQUIRED : exp is an array of integer with dim*nvr allocated space
 *   for the flattened vectors of tropisms. */

int mapcon_coefficients_and_exponents_of_map
 ( int dim, int ind, int nvr, double *cff, int *exp );
/* 
 * DESCRIPTION :
 *   Returns in cff the coefficients of the map of dimension dim
 *   and index ind, where nvr is the number of variables.
 *
 * REQUIRED : cff is an array of doubles with 2*nvr allocated space
 *   for the real and imaginary parts of the coefficients 
 *   and exp is an array of integer with dim*nvr allocated space
 *   for the flattened vectors of tropisms. */

#endif
