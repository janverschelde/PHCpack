/* The prompt_for_setup.h contain the prototype of one function. */

#ifndef __prompt_for_setup_h__
#define __prompt_for_setup_h__

void prompt_for_setup
 ( int *seed, int *dim, int *nbr, int *nva, int *pwr, int *deg, int *vrb,
   int *mode );
/*
 * DESCRIPTION :
 *   Prompts for the parameters to setup an experiment.
 *
 * ON RETURN :
 *   seed     the seed for the random number generator (0 for time);
 *   dim      dimension, total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   nva      number of variables per monomial (0 for random polynomial);
 *   pwr      largest power of a variable;
 *   deg      degree of the series;
 *   vrb      the verbose level;
 *   mode     execution mode, 0, 1 or 2. */

#endif
