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

#endif
