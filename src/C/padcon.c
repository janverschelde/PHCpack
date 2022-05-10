/* The file padcon.c contains the definitions of the functions with
 * prototypes documented in padcon.h. */

#include <stdio.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "padcon.h"

int padcon_set_default_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(735,a,b,c,0);

   return fail;
}

int padcon_clear_parameters ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(736,a,b,c,0);

   return fail;
}

void padcon_write_homotopy_continuation_parameters ( void )
{
   double fltval[2];
   int fail,intval;

   printf("Values of the HOMOTOPY CONTINUATION PARAMETERS :\n");
   fail = padcon_get_homotopy_continuation_parameter(1,fltval);
   printf(" 1. gamma : ");
   printf("%+.14E", fltval[0]); printf("  %+.14E\n", fltval[1]);
   printf(" 2. degree of numerator of Pade approximant    : ");
   fail = padcon_get_homotopy_continuation_parameter(2,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf(" 3. degree of denominator of Pade approximant  : ");
   fail = padcon_get_homotopy_continuation_parameter(3,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf(" 4. maximum step size                          : ");
   fail = padcon_get_homotopy_continuation_parameter(4,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 5. minimum step size                          : ");
   fail = padcon_get_homotopy_continuation_parameter(5,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 6. multiplication factor for the pole radius  : ");
   fail = padcon_get_homotopy_continuation_parameter(6,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 7. multiplication factor for the curvature    : ");
   fail = padcon_get_homotopy_continuation_parameter(7,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 8. tolerance on the residual of the predictor : ");
   fail = padcon_get_homotopy_continuation_parameter(8,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 9. tolerance on the residual of the corrector : ");
   fail = padcon_get_homotopy_continuation_parameter(9,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("10. tolerance on zero series coefficients      : ");
   fail = padcon_get_homotopy_continuation_parameter(10,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("11. maximum number of corrector steps          : ");
   fail = padcon_get_homotopy_continuation_parameter(11,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf("12. maximum steps on a path                    : ");
   fail = padcon_get_homotopy_continuation_parameter(12,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
}

int padcon_write_homotopy_continuation_parameters_to_defined_output ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(874,a,b,c,0);

   return fail;
}

void padcon_tune_homotopy_continuation_parameters ( void )
{
   int choice,intpar,fail;
   double fltpar;
   double gamma[2];

   do
   {
      padcon_write_homotopy_continuation_parameters();
      printf("Type a number to change a value, or 0 to exit : ");
      scanf("%d", &choice);
      if(choice == 1)
      {
         printf("-> give the real part of the new gamma : ");
         scanf("%lf", &gamma[0]);
         printf("-> give the imaginary part of the new gamma : ");
         scanf("%lf", &gamma[1]);
         fail = padcon_set_homotopy_continuation_parameter(1,gamma);
      }
      else if(choice == 2)
      {
         printf("-> give a new numerator degree for the Pade approximant : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 3)
      {
         printf("-> give a new denominator degree for the Pade approximant : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 4)
      {
         printf("-> give a new value for the maximum step size : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 5)
      {
         printf("-> give a new value for the minimum step size  : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 6)
      {
         printf("-> give a new multiplication factor for the pole radius : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 7)
      {
         printf("-> give a new multiplication factor for the curvature : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 8)
      {
         printf("-> give a new tolerance on the predictor residual : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 9)
      {
         printf("-> give a new tolerance on the corrector residual : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 10)
      {
         printf("-> give a new tolerance on a zero series coefficient : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 11)
      {
         printf("-> give a new maximum number of corrector steps : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 12)
      {
         printf("-> give a new maximum number of steps on a path : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      if((choice > 1) && (choice < 13))
      {
         // printf("setting parameter %d to %.3e ...\n", choice, fltpar);
         fail = padcon_set_homotopy_continuation_parameter(choice,&fltpar);
      }
   }
   while (choice != 0);
}

int padcon_get_homotopy_continuation_parameter ( int k, double *val )
{
   int fail,parval;

   if((k == 2) || (k == 3) || (k == 11) || (k == 12))
   {
      fail = _ada_use_c2phc4c(737,&k,&parval,val,0);

      *val = (double) parval; // pass integer value
   }
   else
      fail = _ada_use_c2phc4c(737,&k,&parval,val,0);

   return fail;
}

int padcon_set_homotopy_continuation_parameter ( int k, double *val )
{
   int fail,parval;

   if((k == 2) || (k == 3) || (k == 11) || (k == 12))
   {
      parval = (int) (*val); // pass integer value

      fail = _ada_use_c2phc4c(738,&k,&parval,val,0);
   }
   else
      fail = _ada_use_c2phc4c(738,&k,&parval,val,0);

   return fail;
}

int padcon_reset_homotopy_continuation_parameters ( int prc )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(740,&prc,b,c,0);

   return fail;
}

int padcon_prompt_for_homogenization ( void )
{
   char ans,nlsb;

   printf("\nRunning in homogeneous coordinates ? (y/n) ? ");
   scanf("%c",&ans);
   scanf("%c",&nlsb); /* skip newline symbol */

   if(ans == 'y')
      return 1;
   else
      return 0;
}

int padcon_prompt_for_multi_homogenization ( int nvr )
{
   int choice;

   printf("\n");
   printf("MENU for affine, homogeneous or multi-homogeneous coordinates :\n");
   printf("  0 : in affine coordinates, in the original variables;\n");
   printf("  1 : in 1-homogeous coordinates, in projective space;\n");
   printf("  2 or higher : in multi-homogeous coordinates, in a multi-\n");
   printf("  projective space defined by a partition of the variables.\n");
   printf("Type a number between 0 and %d : ",nvr);

   scanf("%d",&choice);

   return choice;
}

void padcon_define_partition ( int m, int nvr, int *idz )
{
   const int nbc = nvr*10; // 10 characters per symbol
   char smbstr[nbc];       // holds all symbols
   char buffer[10];        // buffer holds the current symbol
   int nbr = nbc;
   int idx,bufidx,setidx;

   int fail = syscon_string_of_symbols(&nbr, smbstr);

   printf("Defining a partition of the variables ...\n");
   bufidx = 0; // index to the buffer
   setidx = 0;
   for(idx=0; idx < nbr; idx++)
   {
      if(smbstr[idx] != ' ')
         buffer[bufidx++] = smbstr[idx];
      else
      {
         buffer[bufidx] = '\0';
         printf("-> which set will %s be in ? ", buffer);
         scanf("%d", &idz[setidx++]);
         bufidx = 0;
      }
   }
   buffer[bufidx] = '\0';
   printf("-> which set will %s be in ? ", buffer);
   scanf("%d", &idz[setidx++]); // last variable

   printf("Your partition :");
   for(idx = 0; idx < nvr; idx++) printf(" %d", idz[idx]);
   printf("\n");
}

int padcon_add_Z0 ( void )
{
   int fail;
   char name[2];

   name[0] = 'Z';
   name[1] = '0';
   fail = syscon_add_symbol(2,name);

   return fail;
}

int padcon_add_symbols ( int m )
{
   char nbr[10];
   char name[20];
   int idx,fail;

   for(idx=1; idx <= m; idx++)
   {
      name[0] = 'Z';
      name[1] = '\0';

      sprintf(nbr,"%d",idx);
      char *ptr = strcat(name,nbr);
      fail = syscon_add_symbol(strlen(name),name);
   }
   return fail;
}

void padcon_standard_projective_transformation ( void )
{
   int fail;

   fail = solcon_standard_one_homogenization();
   fail = copy_container_to_start_solutions();
   fail = copy_target_system_to_container();
   fail = syscon_standard_one_homogenization(0);
   fail = copy_container_to_target_system();
   fail = copy_start_system_to_container();
   fail = syscon_standard_one_homogenization(1);
   fail = copy_container_to_start_system();
   fail = padcon_add_Z0();
}

void padcon_dobldobl_projective_transformation ( void )
{
   int fail;

   fail = solcon_dobldobl_one_homogenization();
   fail = copy_dobldobl_container_to_start_solutions();
   fail = copy_dobldobl_target_system_to_container();
   fail = syscon_dobldobl_one_homogenization(0);
   fail = copy_dobldobl_container_to_target_system();
   fail = copy_dobldobl_start_system_to_container();
   fail = syscon_dobldobl_one_homogenization(1);
   fail = copy_dobldobl_container_to_start_system();
   fail = padcon_add_Z0();
}

void padcon_quaddobl_projective_transformation ( void )
{
   int fail;

   fail = solcon_quaddobl_one_homogenization();
   fail = copy_quaddobl_container_to_start_solutions();
   fail = copy_quaddobl_target_system_to_container();
   fail = syscon_quaddobl_one_homogenization(0);
   fail = copy_quaddobl_container_to_target_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = syscon_quaddobl_one_homogenization(1);
   fail = copy_quaddobl_container_to_start_system();
   fail = padcon_add_Z0();
}

void padcon_standard_multi_projective_transformation
 ( int n, int m, int *idz )
{
   int fail;

   fail = solcon_standard_multi_homogenization(m);
   fail = copy_container_to_start_solutions();
   fail = copy_target_system_to_container();
   fail = syscon_standard_multi_homogenization(n,m,idz,0);
   fail = copy_container_to_target_system();
   fail = copy_start_system_to_container();
   fail = syscon_standard_multi_homogenization(n,m,idz,1);
   fail = copy_container_to_start_system();
   fail = padcon_add_symbols(m);
}

void padcon_dobldobl_multi_projective_transformation
 ( int n, int m, int *idz )
{
   int fail;

   fail = solcon_dobldobl_multi_homogenization(m);
   fail = copy_dobldobl_container_to_start_solutions();
   fail = copy_dobldobl_target_system_to_container();
   fail = syscon_dobldobl_multi_homogenization(n,m,idz,0);
   fail = copy_dobldobl_container_to_target_system();
   fail = copy_dobldobl_start_system_to_container();
   fail = syscon_dobldobl_multi_homogenization(n,m,idz,1);
   fail = copy_dobldobl_container_to_start_system();
   fail = padcon_add_symbols(m);
}

void padcon_quaddobl_multi_projective_transformation
 ( int n, int m, int *idz )
{
   int fail;

   fail = solcon_quaddobl_multi_homogenization(m);
   fail = copy_quaddobl_container_to_start_solutions();
   fail = copy_quaddobl_target_system_to_container();
   fail = syscon_quaddobl_multi_homogenization(n,m,idz,0);
   fail = copy_quaddobl_container_to_target_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = syscon_quaddobl_multi_homogenization(n,m,idz,1);
   fail = copy_quaddobl_container_to_start_system();
   fail = padcon_add_symbols(m);
}

int padcon_standard_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz )
{
   int fail,idx;
   int pars[6];
   int *b;
   double *c;

   pars[0] = 0;   // set precision to double
   pars[1] = nbc;
   pars[2] = verbose;
   pars[3] = mhom;
   pars[4] = locfile;
   pars[5] = nvr;

   if(mhom < 2)
   {
      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,c,0);
      else
      {
         int iname[nbc];
         for(idx=0; idx<nbc; idx++) iname[idx] = (int) name[idx];
         fail = _ada_use_c2phc4c(739,pars,iname,c,0);
      }
   }
   else
   {
      double z[nvr];
      for(int idx=0; idx<nvr; idx++) z[idx] = (double) idz[idx];
 
      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,z,0);
      else
      {
         int iname[nbc];
         for(int i=0; i<nbc; i++) iname[i] = (int) name[i];
         fail = _ada_use_c2phc4c(739,pars,iname,z,0);
      }
   }
   return fail;
}

int padcon_dobldobl_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz )
{
   int fail,idx;
   int pars[6];
   int *b;
   double *c;

   pars[0] = 1;   // set precision to double double
   pars[1] = nbc;
   pars[2] = verbose;
   pars[3] = mhom;
   pars[4] = locfile;
   pars[5] = nvr;

   if(mhom < 2)
   {
      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,c,0);
      else
      {
         int iname[nbc];
         for(idx=0; idx<nbc; idx++) iname[idx] = (int) name[idx];
         fail = _ada_use_c2phc4c(739,pars,iname,c,0);
      }
   }
   else
   {
      double z[nvr];
      for(idx=0; idx<nvr; idx++) z[idx] = (double) idz[idx];

      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,z,0);
      else
      {
         int iname[nbc];
         for(idx=0; idx<nbc; idx++) iname[idx] = (int) name[idx];
         fail = _ada_use_c2phc4c(739,pars,iname,z,0);
      }
   }
   return fail;
}

int padcon_quaddobl_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz )
{
   int fail,idx;
   int pars[6];
   int *b;
   double *c;

   pars[0] = 2;   // set precision to quad double
   pars[1] = nbc;
   pars[2] = verbose;
   pars[3] = mhom;
   pars[4] = locfile;
   pars[5] = nvr;

   if(mhom < 2)
   {
      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,c,0);
      else
      {
         int iname[nbc];
         for(idx=0; idx<nbc; idx++) iname[idx] = (int) name[idx];
         fail = _ada_use_c2phc4c(739,pars,iname,c,0);
      }
   }
   else
   {
      double z[nvr];
      for(idx=0; idx<nvr; idx++) z[idx] = (double) idz[idx];

      if(nbc == 0)
         fail = _ada_use_c2phc4c(739,pars,b,z,0);
      else
      {
         int iname[nbc];
         for(idx=0; idx<nbc; idx++) iname[idx] = (int) name[idx];
         fail = _ada_use_c2phc4c(739,pars,iname,z,0);
      }
   }
   return fail;
}

int padcon_standard_initialize_homotopy ( int verbose, int homo )
{
   int fail;
   int precision = 0;
   int pars[2];
   double *c;

   pars[0] = verbose;
   pars[1] = homo;

   fail = _ada_use_c2phc4c(860,&precision,pars,c,0);

   return fail;
}

int padcon_dobldobl_initialize_homotopy ( int verbose, int homo )
{
   int fail;
   int precision = 1;
   int pars[2];
   double *c;
  
   pars[0] = verbose;
   pars[1] = homo;

   fail = _ada_use_c2phc4c(860,&precision,pars,c,0);

   return fail;
}

int padcon_quaddobl_initialize_homotopy ( int verbose, int homo )
{
   int fail;
   int precision = 2;
   int pars[2];
   double *c;

   pars[0] = verbose;
   pars[1] = homo;

   fail = _ada_use_c2phc4c(860,&precision,pars,c,0);

   return fail;
}

int padcon_standard_initialize_parameter_homotopy ( int idx, int verbose )
{
   int fail;
   int precision = 0;
   int pars[2];
   double *c;

   pars[0] = idx;
   pars[1] = verbose;

   fail = _ada_use_c2phc4c(878,&precision,pars,c,0);

   return fail;
}

int padcon_dobldobl_initialize_parameter_homotopy ( int idx, int verbose )
{
   int fail;
   int precision = 1;
   int pars[2];
   double *c;

   pars[0] = idx;
   pars[1] = verbose;

   fail = _ada_use_c2phc4c(878,&precision,pars,c,0);

   return fail;
}

int padcon_quaddobl_initialize_parameter_homotopy ( int idx, int verbose )
{
   int fail;
   int precision = 2;
   int pars[2];
   double *c;

   pars[0] = idx;
   pars[1] = verbose;

   fail = _ada_use_c2phc4c(878,&precision,pars,c,0);

   return fail;
}

int padcon_initialize_standard_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 0;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c,0);

   return fail;
}

int padcon_initialize_dobldobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 1;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c,0);

   return fail;
}

int padcon_initialize_quaddobl_solution ( int idx, int verbose )
{
   int fail;
   int pars[2];
   double *c;

   pars[0] = 2;    /* double precision */
   pars[1] = idx;

   fail = _ada_use_c2phc4c(861,pars,&verbose,c,0);

   return fail;
}

int padcon_standard_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 0;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c,0);

   *fail = precision;

   return callfail;
}

int padcon_dobldobl_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 1;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c,0);

   *fail = precision;

   return callfail;
}

int padcon_quaddobl_predict_correct ( int* fail, int verbose )
{
   int callfail;
   int precision = 2;
   double *c;

   callfail = _ada_use_c2phc4c(862,&precision,&verbose,c,0);

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

   fail = _ada_use_c2phc4c(863,pars,&verbose,c,0);

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

   fail = _ada_use_c2phc4c(863,pars,&verbose,c,0);

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

   fail = _ada_use_c2phc4c(863,pars,&verbose,c,0);

   if(pars[0] != 0) fail = pars[0];

   return fail;
}

int padcon_get_standard_pole_radius ( double* frp )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp,0);

   return fail;
}

int padcon_get_dobldobl_pole_radius ( double* frp )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp,0);

   return fail;
}

int padcon_get_quaddobl_pole_radius ( double* frp )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp,0);

   return fail;
}

int padcon_get_standard_closest_pole ( double* cre, double* cim )
{
   int fail;
   int precision = 0;
   int *b;
   double nbr[2];

   fail = _ada_use_c2phc4c(866,&precision,b,nbr,0);

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

   fail = _ada_use_c2phc4c(866,&precision,b,nbr,0);

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

   fail = _ada_use_c2phc4c(866,&precision,b,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_standard_t_value ( double *tval )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval,0);

   return fail;
}

int padcon_get_dobldobl_t_value ( double *tval )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval,0);

   return fail;
}

int padcon_get_quaddobl_t_value ( double *tval )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval,0);

   return fail;
}

int padcon_get_standard_step_size ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step,0);

   return fail;
}

int padcon_get_dobldobl_step_size ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step,0);

   return fail;
}

int padcon_get_quaddobl_step_size ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step,0);

   return fail;
}

int padcon_get_standard_series_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step,0);

   return fail;
}

int padcon_get_dobldobl_series_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step,0);

   return fail;
}

int padcon_get_quaddobl_series_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step,0);

   return fail;
}

int padcon_get_standard_pole_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step,0);

   return fail;
}

int padcon_get_dobldobl_pole_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step,0);

   return fail;
}

int padcon_get_quaddobl_pole_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step,0);

   return fail;
}

int padcon_get_standard_estimated_distance ( double *eta )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta,0);

   return fail;
}

int padcon_get_dobldobl_estimated_distance ( double *eta )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta,0);

   return fail;
}

int padcon_get_quaddobl_estimated_distance ( double *eta )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta,0);

   return fail;
}

int padcon_get_standard_hessian_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step,0);

   return fail;
}

int padcon_get_dobldobl_hessian_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step,0);

   return fail;
}

int padcon_get_quaddobl_hessian_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step,0);

   return fail;
}

int padcon_get_standard_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 0;
   inpars[1] = leadidx;
   inpars[2] = idx;

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_dobldobl_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 1;
   inpars[1] = leadidx;
   inpars[2] = idx;

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 2;
   inpars[1] = leadidx;
   inpars[2] = idx;

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_standard_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 0;
   inpars[1] = 1;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_dobldobl_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 1;
   inpars[1] = 1;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 2;
   inpars[1] = 1;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_standard_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 0;
   inpars[1] = 0;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_dobldobl_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 1;
   inpars[1] = 0;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[4];
   double nbr[2];

   inpars[0] = 2;
   inpars[1] = 0;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_standard_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 0;        /* double precision */
   inpars[1] = leadidx;
   inpars[2] = poleidx;

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_dobldobl_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 1;        /* double double precision */
   inpars[1] = leadidx;
   inpars[2] = poleidx;

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr,0);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 2;        /* quad double precision */
   inpars[1] = leadidx;
   inpars[2] = poleidx;

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr,0);

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
   
   fail = _ada_use_c2phc4c(864,&precision,b,c,0);

   return fail;
}

int padcon_clear_dobldobl_data ( void )
{
   int fail;
   int precision = 1;
   int *b;
   double *c;
   
   fail = _ada_use_c2phc4c(864,&precision,b,c,0);

   return fail;
}

int padcon_clear_quaddobl_data ( void )
{
   int fail;
   int precision = 2;
   int *b;
   double *c;
   
   fail = _ada_use_c2phc4c(864,&precision,b,c,0);

   return fail;
}
