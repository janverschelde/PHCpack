/* The file adepath_dd.c contains the definitions of the functions 
 * with prototypes documented in adepath_dd.h. */

#include <stdlib.h>
#include <stdio.h>
#include "adepath_dd.h"
#include "adenewton_dd.h"
#include "adeonepath_dd.h"
#include "ademanypaths_dd.h"

int ade_newton_dd ( int verbose )
{
   int fail = adenewton_dd(verbose);

   return fail;
}

int ade_onepath_dd ( int verbose, double regamma, double imgamma )
{
   int fail = adeonepath_dd(verbose,regamma,imgamma);

   return fail;
}

int ade_manypaths_dd ( int verbose, double regamma, double imgamma )
{
   int fail = ademanypaths_dd(verbose,regamma,imgamma);

   return fail;
}
