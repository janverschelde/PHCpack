/* The file reducers.c contains the definitions of the functions
 * declared in reducers.h. */

#include "reducers.h"

int standard_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(707,&diag,b,c,0);

   return 0;
}

int dobldobl_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(708,&diag,b,c,0);

   return 0;
}

int quaddobl_row_reduce_system ( int diag )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(709,&diag,b,c,0);

   return 0;
}

int standard_nonlinear_reduce_system
 ( int eqmax, int spmax, int rpmax, int *eqcnt, int *spcnt, int *rpcnt )
{
   int fail;
   int max[3];
   int cnt[3];
   double *c;

   max[0] = eqmax;
   max[1] = spmax;
   max[2] = rpmax;

   fail = _ada_use_c2phc4c(710,max,cnt,c,0);

   *eqcnt = cnt[0];
   *spcnt = cnt[1];
   *rpcnt = cnt[2];

   return 0;
}
