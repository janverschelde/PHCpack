/* The file numbtrop.h contains prototypes to the container for managing
 * numerically computed tropisms.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __NUMBTROP_H__
#define __NUMBTROP_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int standard_initialize
 ( int nbt, int dim, int *wnd, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Initializes the numerical tropisms container.
 *
 * ON ENTRY :
 *    nbt    number of tropisms;
 *    dim    length_of_each tropism;
 *    wnd    winding numbers, as many as nbt;
 *    dir    nbt*dim doubles with the coordinates of the tropisms;
 *    err    errors on the tropisms, as many doubles as the value of nbt.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int standard_size ( int *nbt );
/*
 * DESCRIPTION :
 *   Returns in nbt the number of tropisms, stored in standard double
 *   precision, in the numerical tropisms container. */

int standard_retrieve_tropism
 ( int dim, int idx, int *wnd, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Returns one tropisms, stored in standard double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by standard_size.
 *
 * ON RETURN :
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as dim;
 *   err     the error on the tropism. */

int standard_clear( void );
/*
 * DESCRIPTION :
 *   Deallocates the stored numerically computed tropisms,
 *   computed in standard double precision. */

#endif
