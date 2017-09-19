/* The file reducers.c contains the definitions of the functions
 * declared in reducers.h. */

#include "reducers.h"

int standard_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(707,&diag,b,c);

   return 0;
}

int dobldobl_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(708,&diag,b,c);

   return 0;
}

int quaddobl_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(709,&diag,b,c);

   return 0;
}
