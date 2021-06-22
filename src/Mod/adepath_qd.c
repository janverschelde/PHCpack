/* The file adepath_qd.c contains the definitions of the functions 
 * with prototypes documented in adepath_qd.h. */

#include <stdlib.h>
#include <stdio.h>
#include "adepath_qd.h"
#include "adenewton_qd.h"
#include "adeonepath_qd.h"
#include "ademanypaths_qd.h"

int ade_newton_qd ( int verbose )
{
   int fail = adenewton_qd(verbose);

   return fail;
}

int ade_onepath_qd ( int verbose, double regamma, double imgamma )
{
   int fail = adeonepath_qd(verbose,regamma,imgamma);

   return fail;
}

int ade_manypaths_qd ( int verbose, double regamma, double imgamma )
{
   int fail = ademanypaths_qd(verbose,regamma,imgamma);

   return fail;
}
