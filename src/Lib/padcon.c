/* The file padcon.c contains the definitions of the functions with
 * prototypes documented in padcon.h. */

#include "padcon.h"

void padcon_set_default_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(735,a,b,c);
}

void padcon_clear_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(736,a,b,c);
}

int padcon_get_homotopy_continuation_parameter ( int k, double *val )
{
   int fail,parval;

   if((k == 2) || (k == 3) || (k == 11) || (k == 12))
   {
      fail = _ada_use_c2phc4c(737,&k,&parval,val);

      *val = (double) parval; // pass integer value
   }
   else
      fail = _ada_use_c2phc4c(737,&k,&parval,val);

   return fail;
}

int padcon_set_homotopy_continuation_parameter ( int k, double *val )
{
   int fail,parval;

   if((k == 2) || (k == 3) || (k == 11) || (k == 12))
   {
      parval = (int) (*val); // pass integer value

      fail = _ada_use_c2phc4c(738,&k,&parval,val);
   }
   else
      fail = _ada_use_c2phc4c(738,&k,&parval,val);

}
