// The file prompt_newton_setup.h specifies the function to prompt
// for the test parameters on Newton's method.

#ifndef __prompt_newton_setup_h__
#define __prompt_newton_setup_h__

void prompt_newton_setup
 ( int *seed, int *szt, int*nbt, int *dim, int *deg, int *size, int *vrblvl,
   int *mode, int *nbritr, int *nbrcol, int *nbsteps, int *cdata );
/*
 * DESCRIPTION :
 *   Prompts for the parameters to test Newton's method.
 *
 * ON RETURN :
 *   seed      the seed for the random number generator (0 for time);
 *   szt       size of one tile;
 *   nbt       number of tiles, szt*nbt equals the dimension;
 *   dim       the dimension is the number of monomials
 *             and the maximum number of variables in each monomial;
 *   deg       degree of the series;
 *   size      size of the numbers;
 *   vrblvl    verbose level (0 if silent);
 *   mode      execution mode, 0, 1, or 2;
 *   nbritr    number of unimodular multiplications in the making of
 *             the exponent matrix, or type of input system;
 *   nbrcol    number of columns of cyclic n-roots;
 *   nbsteps   the number of Newton steps;
 *   cdata     if 0, then on real data, otherwise on complex data. */

#endif
