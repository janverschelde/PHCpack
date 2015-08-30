/* The file gpupath_dd.c contains the definitions of the functions 
 * with prototypes documented in gpupath_dd.h. */

#include <stdlib.h>
#include <stdio.h>
#include "gpupath_dd.h"
#include "gpunewton_dd.h"
#include "gpuonepath_dd.h"
#include "gpumanypaths_dd.h"

int gpu_newton_dd ( int mode, int verbose )
{
   int fail = gpunewton_dd(mode,verbose);

   return fail;
}

int gpu_onepath_dd ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpuonepath_dd(mode,verbose,regamma,imgamma);

   return fail;
}

int gpu_manypaths_dd ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpumanypaths_dd(mode,verbose,regamma,imgamma);

   return fail;
}
