/* This file "celcon.c" contains the prototypes of the operations
 * declared in "celcon.h". */

#include "celcon.h"

int celcon_read_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(80,a,b,c,0);
   return fail;
}

int celcon_write_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(81,a,b,c,0);
   return fail;
}

int celcon_number_of_cells ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(82,length,b,c,0);
   return fail;
}

int celcon_dimension_of_points ( int *dimension )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(83,dimension,b,c,0);
   return fail;
}

int celcon_is_stable ( int *flag )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(879,flag,b,c,0);

   return fail;
}

int celcon_number_of_original_cells ( int *length )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(880,length,b,c,0);

   return fail;
}

int celcon_number_of_stable_cells ( int *length )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(881,length,b,c,0);

   return fail;
}

int celcon_type_of_mixture ( int *r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(84,r,mix,c,0);
   return fail;
}

int celcon_length_of_supports ( int *r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(85,r,length,c,0);
   return fail;
}

int celcon_get_lifted_point ( int n, int i, int j, double *point )
{
   int fail;
   fail = _ada_use_c2phc4c(86,&i,&j,point,0);
   return fail;
}

int celcon_get_inner_normal ( int n, int i, double *normal )
{
   int *b,fail;
   fail = _ada_use_c2phc4c(87,&i,b,normal,0);
   return fail;
}

int celcon_number_of_points_in_cell ( int i, int r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(88,&i,length,c,0);
   return fail;
}

int celcon_get_point_in_cell ( int n, int i, int j, int k, double *point )
{
   int b[2],fail;
   b[0] = j;
   b[1] = k;
   fail = _ada_use_c2phc4c(89,&i,b,point,0);
   return fail;
}

int celcon_mixed_volume ( int i, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(90,&i,mv,c,0);
   return fail;
}

int celcon_mixed_volume_of_supports ( int *mv )
{
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(239,mv,b,c,0);
   return fail;
}

int celcon_initialize_supports ( int nbr )
{
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(240,&nbr,b,c,0);
   return fail;
}

int celcon_set_type_of_mixture ( int r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(91,&r,mix,c,0);
   return fail;
}

int celcon_append_lifted_point ( int n, int i, double *point )

{
   int fail;
   fail = _ada_use_c2phc4c(92,&i,&n,point,0);
   return fail;
}

int celcon_append_mixed_cell
             ( int n, int r, int k, int labels[], double *normal )
{
   int d[3],fail;
   d[0] = r;
   d[1] = n;
   d[2] = k;
   fail = _ada_use_c2phc4c(93,d,labels,normal,0);
   return fail;
}

int celcon_retrieve_mixed_cell
             ( int n, int r, int i, int labels[], double *normal )
{
   int fail;
   fail = _ada_use_c2phc4c(95,&i,labels,normal,0);
   return fail;
}

int celcon_standard_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(96,a,b,c,0);
   return fail;
}

int celcon_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(460,a,b,c,0);
   return fail;
}

int celcon_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(470,a,b,c,0);
   return fail;
}

int celcon_read_standard_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(97,a,b,c,0);
   return fail;
}

int celcon_read_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(461,a,b,c,0);
   return fail;
}

int celcon_read_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(471,a,b,c,0);
   return fail;
}

int celcon_write_standard_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(98,a,b,c,0);
   return fail;
}

int celcon_write_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(462,a,b,c,0);
   return fail;
}

int celcon_write_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(472,a,b,c,0);
   return fail;
}

int celcon_copy_into_standard_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(99,a,b,c,0);
   return fail;
}

int celcon_copy_into_dobldobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(463,a,b,c,0);
   return fail;
}

int celcon_copy_into_quaddobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(473,a,b,c,0);
   return fail;
}

int celcon_copy_from_standard_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(100,a,b,c,0);
   return fail;
}

int celcon_copy_from_dobldobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(464,a,b,c,0);
   return fail;
}

int celcon_copy_from_quaddobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(474,a,b,c,0);
   return fail;
}

int celcon_standard_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(101,a,b,c,0);
   return fail;
}

int celcon_dobldobl_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(465,a,b,c,0);
   return fail;
}

int celcon_quaddobl_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(475,a,b,c,0);
   return fail;
}

int celcon_solve_standard_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(102,&k,mv,c,0);
   return fail;
}

int celcon_solve_stable_standard_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(882,&k,mv,c,0);
   return fail;
}

int celcon_solve_dobldobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(466,&k,mv,c,0);
   return fail;
}

int celcon_solve_stable_dobldobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(883,&k,mv,c,0);
   return fail;
}

int celcon_solve_quaddobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(476,&k,mv,c,0);
   return fail;
}

int celcon_solve_stable_quaddobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(884,&k,mv,c,0);
   return fail;
}

int celcon_track_standard_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc4c(103,&k,b,c,0);
   return fail;
}

int celcon_track_dobldobl_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc4c(467,&k,b,c,0);
   return fail;
}

int celcon_track_quaddobl_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc4c(477,&k,b,c,0);
   return fail;
}

int celcon_copy_start_standard_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(597,&k,&i,c,0);
   return fail;
}

int celcon_copy_start_dobldobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(598,&k,&i,c,0);
   return fail;
}

int celcon_copy_start_quaddobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(599,&k,&i,c,0);
   return fail;
}

int celcon_copy_target_standard_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(104,&k,&i,c,0);
   return fail;
}

int celcon_copy_target_dobldobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(468,&k,&i,c,0);
   return fail;
}

int celcon_copy_target_quaddobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(478,&k,&i,c,0);
   return fail;
}

int celcon_permute_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(105,a,b,c,0);
   return fail;
}

int celcon_permute_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(469,a,b,c,0);
   return fail;
}

int celcon_permute_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(479,a,b,c,0);
   return fail;
}

int celcon_clear_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(94,a,b,c,0);
   return fail;
}
