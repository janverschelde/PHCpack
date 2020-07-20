/* file product.h contains prototypes for the linear-product
 * root counts and random linear-product systems.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __PRODUCT_H__
#define __PRODUCT_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int supporting_set_structure ( void );
/*
 * DESCRIPTION :
 *   Creates a supporting set structure for the system in the container.
 *
 * REQUIRED : the systems container contains a valid polynomial system. */

int write_set_structure ( void );
/*
 * DESCRIPTION :
 *   Writes the supporting set structure. */

int set_structure_string ( int *nc, char *s );
/*
 * DESCRIPTION :
 *   Returns in s the string representation of the set structure,
 *   with in nc the number of characters in the string. */

int parse_set_structure ( int nc, char *s );
/*
 * DESCRIPTION :
 *   Parses the string s with number of characters nc into a set structure. */

int is_set_structure_supporting ( int *inout );
/*
 * DESCRIPTION :
 *   Verifies whether the stored set structure supports the system in
 *   the container for standard double precision coefficients.
 *   If supporting, then the inout is set to 1 on return, else it is 0.
 *   If on input, inout equals 1, then extra information is written
 *   to screen. */

int linear_product_root_count ( int *r );
/*
 * DESCRIPTION :
 *   Returns in r the root count based on the supporting set structure.
 *
 * REQUIRED : supporting_set_structure() was executed. */

int random_linear_product_system ( void );
/*
 * DESCRIPTION :
 *   Replaces the system in the systems container with 
 *   a random linear-product system based on the supporting set structure.
 *
 * REQUIRED :
 *   supporting_set_structure() was executed and the systems container
 *   still contains the original polynomial system. */

int solve_linear_product_system ( void );
/*
 * DESCRIPTION :
 *   Puts the solution of the random linear-product system in the
 *   solutions container.
 *
 * REQUIRED :
 *   random_linear_product_system() was executed. */

int clear_set_structure ( void );
/*
 * DESCRIPTION :
 *   Clears the set structure constructed with supporting_set_structure. */

int m_homogeneous_Bezout_number ( int *bzn, int *ncp, char *partition );
/*
 * DESCRIPTION :
 *   Returns in bzn a m-homogeneous Bezout number for the system in the
 *   container with standard double precision coefficients. 
 *   Returns in partition the string representation of the partition of
 *   the set of unknowns corresponding the m-homogeneous Bezout number.
 *   The value of ncp on return is the number of characters in partition. */

int m_partition_Bezout_number ( int *bzn, int ncp, char *partition );
/*
 * DESCRIPTION :
 *   Returns in bzn a m-homogeneous Bezout number for the system in the
 *   container with standard double precision coefficients,
 *   for a given partition of the set of unknowns.
 *   The partition is defined in the string partition with ncp characters. */

int m_homogeneous_start_system ( int ncp, char *partition );
/*
 * DESCRIPTION :
 *   Replaces the system in the system container with a random linear-product
 *   start system based on the given partition and the container system.
 *   The number of characters in the string partition equals ncp. */

#endif
