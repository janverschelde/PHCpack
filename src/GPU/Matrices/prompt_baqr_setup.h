/* The prompt_baqr_setup.h contain the prototype of one function. */

#ifndef __prompt_baqr_setup_h__
#define __prompt_baqr_setup_h__

void prompt_baqr_setup
 ( int *seed, int *szt, int *nbt, int *nrows, int *vrb, int *mode );
/*
 * DESCRIPTION :
 *   Prompts for the parameters to test the blocked accelerated
 *   Householder QR.
 *
 * ON RETURN :
 *   seed     the seed for the random number generator (0 for time);
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt equals the number of columns;
 *   nrows    number of rows >= number of columns;
 *   vrb      the verbose level;
 *   mode     execution mode, 0, 1 or 2. */

#endif
