/* The file witsols.c contains the definitions of the functions with
 * prototypes documented in witsols.h. */

#include <stdio.h>
#include "witsols.h"

/* the solvers */

int standard_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(845,pars,&verbose,c);

   return fail;
}

int standard_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(846,pars,&verbose,c);

   return fail;
}

int dobldobl_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(847,pars,&verbose,c);

   return fail;
}

int dobldobl_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(848,pars,&verbose,c);

   return fail;
}

int quaddobl_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(849,pars,&verbose,c);

   return fail;
}

int quaddobl_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose )
{
   int fail = 0;
   int pars[4];
   double *c;

   pars[0] = nbtasks;
   pars[1] = topdim;
   pars[2] = filter;
   pars[3] = factor;

   fail = _ada_use_c2phc4c(850,pars,&verbose,c);

   return fail;
}
