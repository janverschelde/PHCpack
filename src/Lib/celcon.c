/* This file "celcon.c" contains the prototypes of the operations
 * declared in "celcon.h". */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

int celcon_read_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(80,a,b,c);
   return fail;
}

int celcon_write_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(81,a,b,c);
   return fail;
}

int celcon_number_of_cells ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(82,length,b,c);
   return fail;
}

int celcon_dimension_of_points ( int *dimension )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(83,dimension,b,c);
   return fail;
}

int celcon_type_of_mixture ( int *r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(84,r,mix,c);
   return fail;
}

int celcon_length_of_supports ( int *r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(85,r,length,c);
   return fail;
}

int celcon_get_lifted_point ( int n, int i, int j, double *point )
{
   int fail;
   fail = _ada_use_c2phc(86,&i,&j,point);
   return fail;
}

int celcon_get_inner_normal ( int n, int i, double *normal )
{
   int *b,fail;
   fail = _ada_use_c2phc(87,&i,b,normal);
   return fail;
}

int celcon_number_of_points_in_cell ( int i, int r, int *length )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(88,&i,length,c);
   return fail;
}

int celcon_get_point_in_cell ( int n, int i, int j, int k, double *point )
{
   int b[2],fail;
   b[0] = j;
   b[1] = k;
   fail = _ada_use_c2phc(89,&i,b,point);
   return fail;
}

int celcon_mixed_volume ( int i, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(90,&i,mv,c);
   return fail;
}

int celcon_set_type_of_mixture ( int r, int *mix )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(91,&r,mix,c);
   return fail;
}

int celcon_append_lifted_point ( int n, int i, double *point )

{
   int fail;
   fail = _ada_use_c2phc(92,&i,&n,point);
   return fail;
}

int celcon_append_mixed_cell
             ( int n, int r, int k, int labels[k], double *normal )
{
   int d[3],fail;
   d[0] = r;
   d[1] = n;
   d[2] = k;
   fail = _ada_use_c2phc(93,d,labels,normal);
   return fail;
}

int celcon_retrieve_mixed_cell
             ( int n, int r, int i, int labels[], double *normal )
{
   int fail;
   fail = _ada_use_c2phc(95,&i,labels,normal);
   return fail;
}

int celcon_create_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(96,a,b,c);
   return fail;
}

int celcon_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(460,a,b,c);
   return fail;
}

int celcon_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(470,a,b,c);
   return fail;
}

int celcon_read_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(97,a,b,c);
   return fail;
}

int celcon_read_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(461,a,b,c);
   return fail;
}

int celcon_read_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(471,a,b,c);
   return fail;
}

int celcon_write_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(98,a,b,c);
   return fail;
}

int celcon_write_dobldobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(462,a,b,c);
   return fail;
}

int celcon_write_quaddobl_random_coefficient_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(472,a,b,c);
   return fail;
}

int celcon_copy_into_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(99,a,b,c);
   return fail;
}

int celcon_copy_into_dobldobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(463,a,b,c);
   return fail;
}

int celcon_copy_into_quaddobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(473,a,b,c);
   return fail;
}

int celcon_copy_from_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(100,a,b,c);
   return fail;
}

int celcon_copy_from_dobldobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(464,a,b,c);
   return fail;
}

int celcon_copy_from_quaddobl_systems_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(474,a,b,c);
   return fail;
}

int celcon_create_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(101,a,b,c);
   return fail;
}

int celcon_dobldobl_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(465,a,b,c);
   return fail;
}

int celcon_quaddobl_polyhedral_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(475,a,b,c);
   return fail;
}

int celcon_solve_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(102,&k,mv,c);
   return fail;
}

int celcon_solve_dobldobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(466,&k,mv,c);
   return fail;
}

int celcon_solve_quaddobl_start_system ( int k, int *mv )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(476,&k,mv,c);
   return fail;
}

int celcon_track_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc(103,&k,b,c);
   return fail;
}

int celcon_track_dobldobl_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc(467,&k,b,c);
   return fail;
}

int celcon_track_quaddobl_solution_path ( int k, int i, int otp )
{
   int b[2],fail;
   double *c;
   b[0] = i;
   b[1] = otp;
   fail = _ada_use_c2phc(477,&k,b,c);
   return fail;
}

int celcon_copy_target_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(104,&k,&i,c);
   return fail;
}

int celcon_copy_target_dobldobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(468,&k,&i,c);
   return fail;
}

int celcon_copy_target_quaddobl_solution_to_container ( int k, int i )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(478,&k,&i,c);
   return fail;
}

int celcon_permute_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(105,a,b,c);
   return fail;
}

int celcon_permute_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(469,a,b,c);
   return fail;
}

int celcon_permute_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(479,a,b,c);
   return fail;
}

int celcon_clear_mixed_cell_configuration ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(94,a,b,c);
   return fail;
}
