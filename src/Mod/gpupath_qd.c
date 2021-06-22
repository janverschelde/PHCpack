/* The file gpupath_qd.c contains the definitions of the functions 
 * with prototypes documented in gpupath_qd.h. */

#include <stdlib.h>
#include <stdio.h>
#include "gpupath_qd.h"
#include "gpunewton_qd.h"
#include "gpuonepath_qd.h"
#include "gpumanypaths_qd.h"

int gpu_newton_qd ( int mode, int verbose )
{
   int fail = gpunewton_qd(mode,verbose);

   return fail;
}

int gpu_onepath_qd ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpuonepath_qd(mode,verbose,regamma,imgamma);

   return fail;
}

int gpu_manypaths_qd ( int mode, int verbose, double regamma, double imgamma )
{
   int fail = gpumanypaths_qd(mode,verbose,regamma,imgamma);

   return fail;
}
