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

int padcon_standard_track ( int nbc, char* name )
{
   int fail;
   int pars[2];
   int *b;
   double *c;

   pars[0] = 0;   // set precision to double
   pars[1] = nbc;

   if(nbc == 0)
      fail = _ada_use_c2phc4c(739,pars,b,c);
   else
   {
      int iname[nbc];
      for(int i=0; i<nbc; i++) iname[i] = (int) name[i];
      fail = _ada_use_c2phc4c(739,pars,iname,c);
   }
   return fail;
}

int padcon_dobldobl_track ( int nbc, char* name )
{
   int fail;
   int pars[2];
   int *b;
   double *c;

   pars[0] = 1;   // set precision to double double
   pars[1] = nbc;

   if(nbc == 0)
      fail = _ada_use_c2phc4c(739,pars,b,c);
   else
   {
      int iname[nbc];
      for(int i=0; i<nbc; i++) iname[i] = (int) name[i];
      fail = _ada_use_c2phc4c(739,pars,iname,c);
   }
   return fail;
}

int padcon_quaddobl_track ( int nbc, char* name )
{
   int fail;
   int pars[2];
   int *b;
   double *c;

   pars[0] = 2;   // set precision to quad double
   pars[1] = nbc;

   if(nbc == 0)
      fail = _ada_use_c2phc4c(739,pars,b,c);
   else
   {
      int iname[nbc];
      for(int i=0; i<nbc; i++) iname[i] = (int) name[i];
      fail = _ada_use_c2phc4c(739,pars,iname,c);
   }
   return fail;
}
