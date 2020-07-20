/* The file sweep.c contains the definitions of the functions
 * declared in sweep.h. */

#include "sweep.h"

int sweep_define_parameters_numerically ( int nq, int nv, int np, int *pars )
{
   int fail = 0;
   int nqvp[3];
   double *c;

   nqvp[0] = nq;
   nqvp[1] = nv;
   nqvp[2] = np;
   fail = _ada_use_c2phc4c(610,nqvp,pars,c,0);

   return fail;
}

int sweep_define_parameters_symbolically
 ( int nq, int nv, int np, int nc, char *pars )
{
   int k,fail = 0;
   int nqvp[4];
   int idx[np];
   double *c;

   nqvp[0] = nq;
   nqvp[1] = nv;
   nqvp[2] = np;
   nqvp[3] = nc;
   for(k=0; k<np; k++) idx[k] = (int) pars[k];

   fail = _ada_use_c2phc4c(611,nqvp,idx,c,0);

   return fail;
}

int sweep_get_number_of_equations ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(612,n,b,c,0);

   return fail;
}

int sweep_get_number_of_variables ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(613,n,b,c,0);

   return fail;
}

int sweep_get_number_of_parameters ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(614,n,b,c,0);

   return fail;
}

int sweep_get_indices_numerically ( int *idx )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(615,idx,b,c,0);

   return fail;
}

int sweep_get_indices_symbolically ( int *nc, char *pars )
{
   int fail,np;
   double *c;

   fail = sweep_get_number_of_parameters(&np);
   {
      int k,buf[np*20]; /* assumes no more than 20 characters per symbol */

      fail = _ada_use_c2phc4c(617,nc,buf,c,0);

      for(k=0; k<(*nc); k++) pars[k] = (char) buf[k];
   }
   return fail;
}

int sweep_clear_definitions ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(617,a,b,c,0);

   return fail;
}

int sweep_set_standard_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 0; /* standard precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_set_standard_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 0; /* standard precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_set_dobldobl_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 1; /* double double precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_set_dobldobl_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 1; /* double double precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_set_quaddobl_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 2; /* quad double precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_set_quaddobl_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 2; /* quad double precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(618,ptos,&n,c,0);

   return fail;
}

int sweep_get_standard_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 0; /* standard precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_get_standard_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 0; /* standard precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_get_dobldobl_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 1; /* double double precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_get_dobldobl_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 1; /* double double precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_get_quaddobl_start ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 2; /* quad double precision */
   ptos[1] = 0; /* start values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_get_quaddobl_target ( int n, double *c )
{
   int fail,ptos[2];

   ptos[0] = 2; /* quad double precision */
   ptos[1] = 1; /* target values */

   fail = _ada_use_c2phc4c(619,ptos,&n,c,0);

   return fail;
}

int sweep_standard_complex_run
 ( int gchoice, double *regamma, double *imgamma )
{
   int fail;
   int pg[2];
   int *b;
   double c[2];

   pg[0] = 0; /* standard double precision */
   pg[1] = gchoice; /* type of choice for gamma */
   if(gchoice == 2)
   {
      c[0] = *regamma;
      c[1] = *imgamma;
   }
   fail = _ada_use_c2phc4c(620,pg,b,c,0);

   return fail;
}

int sweep_dobldobl_complex_run
 ( int gchoice, double *regamma, double *imgamma )
{
   int fail;
   int pg[2];
   int *b;
   double c[4];

   pg[0] = 1; /* double double precision */
   pg[1] = gchoice;   /* type of choice for gamma */
   if(gchoice == 2)
   {
      c[0] = regamma[0]; c[1] = regamma[1]; /* regamma is a double double */
      c[2] = imgamma[0]; c[3] = imgamma[1]; /* imgamma is a double double */
   }
   fail = _ada_use_c2phc4c(620,pg,b,c,0);

   return fail;
}

int sweep_quaddobl_complex_run
 ( int gchoice, double *regamma, double *imgamma )
{
   int fail;
   int pg[2];
   int *b;
   double c[8];

   pg[0] = 2; /* quad double precision */
   pg[1] = gchoice;   /* type of choice for gamma */
   if(gchoice == 2)
   {
      c[0] = regamma[0]; c[1] = regamma[1]; /* regamma is a quad double */
      c[2] = regamma[2]; c[3] = regamma[3];
      c[0] = imgamma[0]; c[1] = imgamma[1]; /* imgamma is a quad double */
      c[2] = imgamma[2]; c[3] = imgamma[3];
   }
   fail = _ada_use_c2phc4c(620,pg,b,c,0);

   return fail;
}

int sweep_standard_real_run ( void )
{
   int fail,*b;
   double *c;
   int precision = 0;

   fail = _ada_use_c2phc4c(621,&precision,b,c,0);

   return fail;
}

int sweep_dobldobl_real_run ( void )
{
   int fail,*b;
   double *c;
   int precision = 1;

   fail = _ada_use_c2phc4c(621,&precision,b,c,0);

   return fail;
}

int sweep_quaddobl_real_run ( void )
{
   int fail,*b;
   double *c;
   int precision = 2;

   fail = _ada_use_c2phc4c(621,&precision,b,c,0);

   return fail;
}
