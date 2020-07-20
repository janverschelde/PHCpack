/* The file multiplicity.c contains the definitions of the functions
 * declared in multiplicity.h. */

#include "multiplicity.h"

int standard_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert )
{
   int fail = 0;

   *multiplicity = order;
   hilbert[0] = verbose;

   fail = _ada_use_c2phc4c(732,multiplicity,hilbert,&tol,0);

   return fail;
}

int dobldobl_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert )
{
   int fail = 0;

   *multiplicity = order;
   hilbert[0] = verbose;

   fail = _ada_use_c2phc4c(733,multiplicity,hilbert,&tol,0);

   return fail;
}

int quaddobl_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert )
{
   int fail = 0;

   *multiplicity = order;
   hilbert[0] = verbose;

   fail = _ada_use_c2phc4c(734,multiplicity,hilbert,&tol,0);

   return fail;
}
