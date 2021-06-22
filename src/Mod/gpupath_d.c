/* The file gpupath_d.c contains the definitions of the functions 
 * with prototypes documented in gpupath_d.h. */

#include <stdlib.h>
#include <stdio.h>
#include "gpupath_d.h"
#include "gpunewton_d.h"
#include "gpuonepath_d.h"
#include "gpumanypaths_d.h"

int gpu_newton_d ( int mode, int verbose )
{
   int fail = gpunewton_d(mode,verbose);

   return fail;
}

int gpu_onepath_d ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpuonepath_d(mode,verbose,regamma,imgamma);

   return fail;
}

int gpu_manypaths_d ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpumanypaths_d(mode,verbose,regamma,imgamma);

   return fail;
}
