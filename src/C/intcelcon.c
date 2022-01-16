/* This file "intcelcon.c" contains the prototypes of the operations
 * declared in "intcelcon.h". */

#include "intcelcon.h"

int intcelcon_read_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(741,a,b,c,0);
   return fail;
}

int intcelcon_write_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(742,a,b,c,0);
   return fail;
}

int intcelcon_number_of_cells ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(743,length,b,c,0);
   return fail;
}

int intcelcon_dimension_of_points ( int *dimension )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(744,dimension,b,c,0);
   return fail;
}

int intcelcon_type_of_mixture ( int *r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(745,r,mix,c,0);
   return fail;
}

int intcelcon_length_of_supports ( int *r, int *length )
{
   int fail,i;
   double *c;
   fail = _ada_use_c2phc4c(746,r,length,c,0);
   return fail;
}

int intcelcon_get_lifted_point ( int n, int i, int j, int *point )
{
   int fail,k;
   double fltpoint[n];

   fail = _ada_use_c2phc4c(747,&i,&j,fltpoint,0);

   for(k=0; k<n; k++) point[k] = (int) fltpoint[k];

   return fail;
}

int intcelcon_get_inner_normal ( int n, int i, int *normal )
{
   int *b,fail,k;
   double dblnormal[n];

   fail = _ada_use_c2phc4c(748,&i,b,dblnormal,0);

   for(k=0; k<n; k++) normal[k] = (int) dblnormal[k];

   return fail;
}

int intcelcon_number_of_points_in_cell ( int i, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(749,&i,length,c,0);
   return fail;
}

int intcelcon_get_point_in_cell ( int n, int i, int j, int k, int *point )
{
   int b[2],fail,kk;
   double dblpoint[n];

   b[0] = j;
   b[1] = k;
   fail = _ada_use_c2phc4c(750,&i,b,dblpoint,0);

   for(kk=0; kk<n; kk++) point[kk] = (int) dblpoint[kk];

   return fail;
}

int intcelcon_mixed_volume ( int i, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(751,&i,mv,c,0);
   return fail;
}

int intcelcon_initialize_supports ( int nbr )
{
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(752,&nbr,b,c,0);
   return fail;
}

int intcelcon_set_type_of_mixture ( int r, int *mix )
{
   int fail,i;
   double *c;
   fail = _ada_use_c2phc4c(753,&r,mix,c,0);
   return fail;
}

int intcelcon_append_lifted_point ( int n, int i, int *point )
{
   int fail,k;
   double dblpoint[n];

   for(k=0; k<n; k++) dblpoint[k] = (double) point[k];

   fail = _ada_use_c2phc4c(754,&i,&n,dblpoint,0);

   return fail;
}

int intcelcon_append_mixed_cell
 ( int n, int r, int k, int labels[], int *normal )
{
   int d[3],fail,i;
   double dblnormal[n];

   for(i=0; i<n; i++) dblnormal[i] = (double) normal[i];

   d[0] = r;
   d[1] = n;
   d[2] = k;
   fail = _ada_use_c2phc4c(754,d,labels,dblnormal,0);

   return fail;
}

int intcelcon_retrieve_mixed_cell
 ( int n, int r, int i, int labels[], int *normal )
{
   int fail,k;
   double dblnormal[n];

   fail = _ada_use_c2phc4c(756,&i,labels,dblnormal,0);

   for(k=0; k<n; k++) normal[k] = (int) dblnormal[k];

   return fail;
}

int intcelcon_make_subdivision ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(758,a,b,c,0);

   return fail;
}

int intcelcon_clear_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(755,a,b,c,0);
   return fail;
}
