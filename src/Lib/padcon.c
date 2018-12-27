/* The file padcon.c contains the definitions of the functions with
 * prototypes documented in padcon.h. */

#include "padcon.h"

int padcon_set_default_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(735,a,b,c);

   return fail;
}

int padcon_clear_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(736,a,b,c);

   return fail;
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

int padcon_standard_track ( int nbc, char* name, int verbose )
{
   int fail;
   int pars[3];
   int *b;
   double *c;

   pars[0] = 0;   // set precision to double
   pars[1] = nbc;
   pars[2] = verbose;

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

int padcon_dobldobl_track ( int nbc, char* name, int verbose )
{
   int fail;
   int pars[3];
   int *b;
   double *c;

   pars[0] = 1;   // set precision to double double
   pars[1] = nbc;
   pars[2] = verbose;

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

int padcon_quaddobl_track ( int nbc, char* name, int verbose )
{
   int fail;
   int pars[3];
   int *b;
   double *c;

   pars[0] = 2;   // set precision to quad double
   pars[1] = nbc;
   pars[2] = verbose;

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

int padcon_standard_initialize_homotopy ( int verbose )
{
   int fail;
   int precision = 0;
   double *c;

   fail = _ada_use_c2phc4c(860,&precision,&verbose,c);

   return fail;
}

int padcon_dobldobl_initialize_homotopy ( int verbose )
{
   int fail;
   int precision = 1;
   double *c;

   fail = _ada_use_c2phc4c(860,&precision,&verbose,c);

   return fail;
}

int padcon_quaddobl_initialize_homotopy ( int verbose )
{
   int fail;
   int precision = 2;
   double *c;

   fail = _ada_use_c2phc4c(860,&precision,&verbose,c);

   return fail;
}

int padcon_initialize_standard_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 0;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c);

   return fail;
}

int padcon_initialize_dobldobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 1;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c);

   return fail;
}

int padcon_initialize_quaddobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 2;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c);

   return fail;
}

int padcon_standard_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 0;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c);

   *fail = precision;

   return callfail;
}

int padcon_dobldobl_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 1;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c);

   *fail = precision;

   return callfail;
}

int padcon_quaddobl_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 2;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c);

   *fail = precision;

   return callfail;
}

int padcon_get_standard_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 0;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(863,pars,&verbose,c);

   if(pars[0] != 0) fail = pars[0];

   return fail;
}

int padcon_get_dobldobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 1;    /* double double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(863,pars,&verbose,c);

   if(pars[0] != 0) fail = pars[0];

   return fail;
}

int padcon_get_quaddobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 2;    /* quad double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(863,pars,&verbose,c);

   if(pars[0] != 0) fail = pars[0];

   return fail;
}

int padcon_get_standard_forward_pole_radius ( double* frp )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp);

   return fail;
}

int padcon_get_dobldobl_forward_pole_radius ( double* frp )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp);

   return fail;
}

int padcon_get_quaddobl_forward_pole_radius ( double* frp )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp);

   return fail;
}

int padcon_get_standard_closest_pole ( double* cre, double* cim )
{
   int fail;
   int precision = 0;
   int *b;
   double nbr[2];

   fail = _ada_use_c2phc4c(866,&precision,b,nbr);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_dobldobl_closest_pole ( double* cre, double* cim )
{
   int fail;
   int precision = 1;
   int *b;
   double nbr[2];

   fail = _ada_use_c2phc4c(866,&precision,b,nbr);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_closest_pole ( double* cre, double* cim )
{
   int fail;
   int precision = 2;
   int *b;
   double nbr[2];

   fail = _ada_use_c2phc4c(866,&precision,b,nbr);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_clear_standard_data ( void )
{
   int fail;
   int precision = 0;
   int *b;
   double *c;
   
   fail = _ada_use_c2phc4c(864,&precision,b,c);

   return fail;
}

int padcon_clear_dobldobl_data ( void )
{
   int fail;
   int precision = 1;
   int *b;
   double *c;
   
   fail = _ada_use_c2phc4c(864,&precision,b,c);

   return fail;
}

int padcon_clear_quaddobl_data ( void )
{
   int fail;
   int precision = 2;
   int *b;
   double *c;
   
   fail = _ada_use_c2phc4c(864,&precision,b,c);

   return fail;
}
