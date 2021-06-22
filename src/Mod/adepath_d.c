/* The file adepath_d.c contains the definitions of the functions 
 * with prototypes documented in adepath_d.h. */

#include <stdlib.h>
#include <stdio.h>
#include "adepath_d.h"
#include "adenewton_d.h"
#include "adeonepath_d.h"
#include "ademanypaths_d.h"

int ade_newton_d ( int verbose )
{
   int fail = adenewton_d(verbose);

   return fail;
}

int ade_onepath_d ( int verbose, double regamma, double imgamma )
{
   int fail = adeonepath_d(verbose,regamma,imgamma);

   return fail;
}

int ade_manypaths_d ( int verbose, double regamma, double imgamma )
{
   int fail = ademanypaths_d(verbose,regamma,imgamma);

   return fail;
}
