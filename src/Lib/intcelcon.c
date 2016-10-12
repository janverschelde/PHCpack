/* This file "intcelcon.c" contains the prototypes of the operations
 * declared in "intcelcon.h". */

#include "intcelcon.h"

int intcelcon_read_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(741,a,b,c);
   return fail;
}

int intcelcon_write_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(742,a,b,c);
   return fail;
}

int intcelcon_number_of_cells ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(743,length,b,c);
   return fail;
}

int intcelcon_dimension_of_points ( int *dimension )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(744,dimension,b,c);
   return fail;
}

int intcelcon_type_of_mixture ( int *r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(745,r,mix,c);
   return fail;
}

int intcelcon_length_of_supports ( int *r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(746,r,length,c);
   return fail;
}

int intcelcon_get_lifted_point ( int n, int i, int j, int *point )
{
   int fail,k;
   double fltpoint[n];

   fail = _ada_use_c2phc(747,&i,&j,fltpoint);

   for(k=0; k<n; i++) point[k] = (int) fltpoint[k];

   return fail;
}

int intcelcon_get_inner_normal ( int n, int i, int *normal )
{
   int *b,fail,k;
   double dblnormal[n];

   fail = _ada_use_c2phc(748,&i,b,dblnormal);

   for(k=0; k<n; k++) normal[k] = (int) dblnormal[k];

   return fail;
}

int intcelcon_number_of_points_in_cell ( int i, int r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(749,&i,length,c);
   return fail;
}

int intcelcon_get_point_in_cell ( int n, int i, int j, int k, int *point )
{
   int b[2],fail,kk;
   double dblpoint[n];

   b[0] = j;
   b[1] = k;
   fail = _ada_use_c2phc(750,&i,b,dblpoint);

   for(kk=0; kk<n; kk++) point[kk] = (int) dblpoint[kk];

   return fail;
}

int intcelcon_mixed_volume ( int i, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(751,&i,mv,c);
   return fail;
}

int intcelcon_initialize_supports ( int nbr )
{
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc(752,&nbr,b,c);
   return fail;
}

int intcelcon_set_type_of_mixture ( int r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(753,&r,mix,c);
   return fail;
}

int intcelcon_append_lifted_point ( int n, int i, int *point )
{
   int fail,k;
   double dblpoint[n];

   for(k=0; k<n; k++) dblpoint[k] = (double) point[k];

   fail = _ada_use_c2phc(754,&i,&n,dblpoint);

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
   fail = _ada_use_c2phc(754,d,labels,dblnormal);

   return fail;
}

int intcelcon_retrieve_mixed_cell
 ( int n, int r, int i, int labels[], int *normal )
{
   int fail,k;
   double dblnormal[n];

   fail = _ada_use_c2phc(756,&i,labels,dblnormal);

   for(k=0; k<n; k++) normal[k] = (int) dblnormal[k];

   return fail;
}

int intcelcon_clear_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(755,a,b,c);
   return fail;
}
