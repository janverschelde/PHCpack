/* The file numbtrop.c contains the definitions of the functions with
 * prototypes documented in numbtrop.h. */

#include "numbtrop.h"

int numbtrop_standard_initialize
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;
   for(k=0; k<nbt*dim; k++) cff[k] = dir[k];
   for(k=0; k<nbt; k++) cff[nbt*dim+k] = err[k];

   fail = _ada_use_c2phc4c(711,nd,wnd,cff,0);

   return fail;
}

int numbtrop_dobldobl_initialize
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[2*nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;
   for(k=0; k<2*nbt*dim; k++) cff[k] = dir[k];
   for(k=0; k<2*nbt; k++) cff[2*nbt*dim+k] = err[k];

   fail = _ada_use_c2phc4c(712,nd,wnd,cff,0);

   return fail;
}

int numbtrop_quaddobl_initialize
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[4*nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;
   for(k=0; k<4*nbt*dim; k++) cff[k] = dir[k];
   for(k=0; k<4*nbt; k++) cff[4*nbt*dim+k] = err[k];

   fail = _ada_use_c2phc4c(713,nd,wnd,cff,0);

   return fail;
}

int numbtrop_standard_retrieve
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;

   fail = _ada_use_c2phc4c(717,nd,wnd,cff,0);

   for(k=0; k<nbt*dim; k++) dir[k] = cff[k];
   for(k=0; k<nbt; k++) err[k] = cff[nbt*dim+k];

   return fail;
}

int numbtrop_dobldobl_retrieve
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[2*nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;

   fail = _ada_use_c2phc4c(718,nd,wnd,cff,0);

   for(k=0; k<2*nbt*dim; k++) dir[k] = cff[k];
   for(k=0; k<2*nbt; k++) err[k] = cff[2*nbt*dim+k];

   return fail;
}

int numbtrop_quaddobl_retrieve
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[4*nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;

   fail = _ada_use_c2phc4c(719,nd,wnd,cff,0);

   for(k=0; k<4*nbt*dim; k++) dir[k] = cff[k];
   for(k=0; k<4*nbt; k++) err[k] = cff[4*nbt*dim+k];

   return fail;
}

int numbtrop_standard_size ( int *nbt )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(720,nbt,b,c,0);

   return fail;
}

int numbtrop_dobldobl_size ( int *nbt )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(721,nbt,b,c,0);

   return fail;
}

int numbtrop_quaddobl_size ( int *nbt )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(722,nbt,b,c,0);

   return fail;
}

int numbtrop_standard_dimension ( int *dim )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(729,dim,b,c,0);

   return fail;
}

int numbtrop_dobldobl_dimension ( int *dim )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(730,dim,b,c,0);

   return fail;
}

int numbtrop_quaddobl_dimension ( int *dim )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(731,dim,b,c,0);

   return fail;
}

int numbtrop_store_standard_tropism
 ( int dim, int idx, int wnd, double *dir, double err )
{
   int fail,k;
   int nd[2];
   double cff[dim+1];

   nd[0] = dim; nd[1] = idx;
   for(k=0; k<dim; k++) cff[k] = dir[k];
   cff[dim] = err;

   fail = _ada_use_c2phc4c(714,nd,&wnd,cff,0);

   return fail;
}

int numbtrop_store_dobldobl_tropism
 ( int dim, int idx, int wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[2*(dim+1)];

   nd[0] = dim; nd[1] = idx;
   for(k=0; k<2*dim; k++) cff[k] = dir[k];
   cff[2*dim] = err[0];
   cff[2*dim+1] = err[1];

   fail = _ada_use_c2phc4c(715,nd,&wnd,cff,0);

   return fail;
}

int numbtrop_store_quaddobl_tropism
 ( int dim, int idx, int wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[4*(dim+1)];

   nd[0] = dim; nd[1] = idx;
   for(k=0; k<4*dim; k++) cff[k] = dir[k];
   cff[4*dim] = err[0];
   cff[4*dim+1] = err[1];
   cff[4*dim+2] = err[2];
   cff[4*dim+3] = err[3];

   fail = _ada_use_c2phc4c(716,nd,&wnd,cff,0);

   return fail;
}

int numbtrop_standard_retrieve_tropism
 ( int dim, int idx, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[dim+1];

   nd[0] = dim;
   nd[1] = idx;

   fail = _ada_use_c2phc4c(723,nd,wnd,cff,0);

   for(k=0; k<dim; k++) dir[k] = cff[k];
   *err = cff[dim];

   return fail;
}

int numbtrop_dobldobl_retrieve_tropism
 ( int dim, int idx, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[2*(dim+1)];

   nd[0] = dim;
   nd[1] = idx;

   fail = _ada_use_c2phc4c(724,nd,wnd,cff,0);

   for(k=0; k<2*dim; k++) dir[k] = cff[k];
   err[0] = cff[2*dim];
   err[1] = cff[2*dim+1];

   return fail;
}

int numbtrop_quaddobl_retrieve_tropism
 ( int dim, int idx, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[4*(dim+1)];

   nd[0] = dim;
   nd[1] = idx;

   fail = _ada_use_c2phc4c(725,nd,wnd,cff,0);

   for(k=0; k<4*dim; k++) dir[k] = cff[k];
   err[0] = cff[4*dim];
   err[1] = cff[4*dim+1];
   err[2] = cff[4*dim+2];
   err[3] = cff[4*dim+3];

   return fail;
}

int numbtrop_standard_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(726,a,b,c,0);

   return fail;
}

int numbtrop_dobldobl_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(727,a,b,c,0);

   return fail;
}

int numbtrop_quaddobl_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(728,a,b,c,0);

   return fail;
}
