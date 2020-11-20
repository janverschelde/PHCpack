/* The file giftwrappers.c contains the definitions of the functions 
 * with prototypes documented in giftwrappers.h. */

#include <stdlib.h>
#include <stdio.h>
#include "giftwrappers.h"

int convex_hull_2d ( int nc_pts, char *pts, int *nc_hull, char *hull )
{
   int fail,i;
   int points[10*nc_pts];
   double *c;

   for(i=0; i<nc_pts; i++) points[i] = (int) pts[i];

   fail = _ada_use_c2phc4c(580,&nc_pts,points,c,0);

   *nc_hull = nc_pts;

   for(i=0; i<nc_pts; i++) hull[i] = (char) points[i];
  
   return fail;
}

int convex_hull ( int nc_pts, char *pts )
{
   int fail,i;
   int points[nc_pts];
   double *c;

   for(i=0; i<nc_pts; i++) points[i] = (int) pts[i];

   fail = _ada_use_c2phc4c(581,&nc_pts,points,c,0);

   return fail;
}

int number_of_facets ( int dim, int *nbr )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(582,&dim,nbr,c,0);

   return fail;
}

int retrieve_facet ( int dim, int fcn, int *nc_rep, char *fctrep )
{
   int fail,k;
   int a = dim;
   int b[100*dim];
   double *c;

   b[0] = fcn;

   fail = _ada_use_c2phc4c(583,&a,b,c,0);

   *nc_rep = a;
   for(k=0; k<a; k++) fctrep[k] = (char) b[k];
   fctrep[a] = '\0';

   return fail;
}

int clear_3d_facets ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(584,a,b,c,0);

   return fail;
}

int clear_4d_facets ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(585,a,b,c,0);

   return fail;
}

int support_size ( int idx ) 
{
   int fail,nbr;
   int *b;
   double *c;

   nbr = idx;

   fail = _ada_use_c2phc4c(586,&nbr,b,c,0);

   return nbr;
}

int support_string ( int size, char *supp )
{
   int fail,i;
   int support[size];
   double *c;

   fail = _ada_use_c2phc4c(587,&size,support,c,0);
   for(i=0; i<size; i++) supp[i] = (char) support[i];
   supp[size] = '\0';

   return fail;
}

int clear_support_string ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(588,a,b,c,0);

   return fail;
}

int initial_form ( int dim, int nbc, char *normal )
{
   int fail,i;
   int pars[2];
   int nrm[nbc];
   double *c;

   pars[0] = dim;
   pars[1] = nbc;
   /* printf("the normal in C is %s\n", normal); */
   for(i=0; i<nbc; i++) nrm[i] = (int) normal[i];
   fail = _ada_use_c2phc4c(589,pars,nrm,c,0);

   return fail;
}
