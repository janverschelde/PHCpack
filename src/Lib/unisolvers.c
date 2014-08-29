/* file unisolvers.c contains the definitions of the functions 
 * with prototypes documented in unisolvers.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern void adainit ( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal ( void );

int solve_with_standard_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc(272,&max,nit,&eps);

   return fail;
}

int solve_with_double_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc(273,&max,nit,&eps);

   return fail;
}

int solve_with_quad_doubles ( int max, double eps, int *nit )
{
   int fail = _ada_use_c2phc(274,&max,nit,&eps);

   return fail;
}
