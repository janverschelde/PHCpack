/* file unisolvers.c contains the definitions of the functions 
 * with prototypes documented in unisolvers.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "unisolvers.h"

int solve_with_standard_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc4c(272,&max,nit,&eps,0);

   return fail;
}

int solve_with_double_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc4c(273,&max,nit,&eps,0);

   return fail;
}

int solve_with_quad_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc4c(274,&max,nit,&eps,0);

   return fail;
}

int solve_with_multiprecision ( int dcp, int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc4c(275,&max,nit,&eps,0);

   return fail;
}
