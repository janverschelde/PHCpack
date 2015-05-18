/* This file "mapcon.c" contains the definitions of the operations
 * declared in the file "mapcon.h". */

#include "mapcon.h"

int mapcon_solve_system ( int puretopdim )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(430,&puretopdim,b,c);
   return fail;
}

int mapcon_write_maps ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(431,a,b,c);
   return fail;
}

int mapcon_clear_maps ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(432,a,b,c);
   return fail;
}

int mapcon_top_dimension ( int *dim )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(433,dim,b,c);
   return fail;
}

int mapcon_number_of_maps ( int dim, int *nbmaps )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(434,&dim,nbmaps,c);
   return fail;
}

int mapcon_degree_of_map ( int dim, int ind, int *deg )
{
   int a[2], fail;
   double *c;

   a[0] = dim;
   a[1] = ind;
   fail = _ada_use_c2phc(435,a,deg,c);

   return fail;
}

int mapcon_coefficients_of_map ( int dim, int ind, int nvr, double *cff )
{
   int a[3], *b, fail;
  
   a[0] = dim;
   a[1] = ind;
   a[2] = nvr;
   fail = _ada_use_c2phc(436,a,b,cff);

   return fail;
}

int mapcon_exponents_of_map ( int dim, int ind, int nvr, int *exp )
{
   int a[3], fail;
   double *c;
  
   a[0] = dim;
   a[1] = ind;
   a[2] = nvr;
   fail = _ada_use_c2phc(437,a,exp,c);

   return fail;
}

int mapcon_coefficients_and_exponents_of_map
 ( int dim, int ind, int nvr, double *cff, int *exp )
{
   int a[3], fail;
  
   a[0] = dim;
   a[1] = ind;
   a[2] = nvr;
   fail = _ada_use_c2phc(438,a,exp,cff);

   return fail;
}
