/* The file padcon.c contains the definitions of the functions with
 * prototypes documented in padcon.h. */

#include <stdio.h>
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
   printf(" 6. multiplication factor for the series step  : ");
   fail = padcon_get_homotopy_continuation_parameter(6,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 7. multiplication factor for the pole radius  : ");
   fail = padcon_get_homotopy_continuation_parameter(7,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 8. multiplication factor for the curvature    : ");
   fail = padcon_get_homotopy_continuation_parameter(8,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 9. tolerance on the residual of the predictor : ");
   fail = padcon_get_homotopy_continuation_parameter(9,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("10. tolerance on the residual of the corrector : ");
   fail = padcon_get_homotopy_continuation_parameter(10,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("11. tolerance on zero series coefficients      : ");
   fail = padcon_get_homotopy_continuation_parameter(11,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("12. maximum number of corrector steps          : ");
   fail = padcon_get_homotopy_continuation_parameter(12,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf("13. maximum steps on a path                    : ");
   fail = padcon_get_homotopy_continuation_parameter(13,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
}

int padcon_write_homotopy_continuation_parameters_to_defined_output ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(874,a,b,c);

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

   return fail;
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

int padcon_standard_initialize_parameter_homotopy ( int idx, int verbose )
{
   int fail;
   int precision = 0;
   int pars[2];
   double *c;

   pars[0] = idx;
   pars[1] = verbose;

   fail = _ada_use_c2phc4c(878,&precision,pars,c);

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

   fail = _ada_use_c2phc4c(878,&precision,pars,c);

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

   fail = _ada_use_c2phc4c(878,&precision,pars,c);

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

int padcon_get_standard_pole_radius ( double* frp )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp);

   return fail;
}

int padcon_get_dobldobl_pole_radius ( double* frp )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(865,&precision,b,frp);

   return fail;
}

int padcon_get_quaddobl_pole_radius ( double* frp )
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

int padcon_get_standard_t_value ( double *tval )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval);

   return fail;
}

int padcon_get_dobldobl_t_value ( double *tval )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval);

   return fail;
}

int padcon_get_quaddobl_t_value ( double *tval )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(867,&precision,b,tval);

   return fail;
}

int padcon_get_standard_step_size ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step);

   return fail;
}

int padcon_get_dobldobl_step_size ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step);

   return fail;
}

int padcon_get_quaddobl_step_size ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(868,&precision,b,step);

   return fail;
}

int padcon_get_standard_series_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step);

   return fail;
}

int padcon_get_dobldobl_series_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step);

   return fail;
}

int padcon_get_quaddobl_series_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(885,&precision,b,step);

   return fail;
}

int padcon_get_standard_pole_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step);

   return fail;
}

int padcon_get_dobldobl_pole_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step);

   return fail;
}

int padcon_get_quaddobl_pole_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(886,&precision,b,step);

   return fail;
}

int padcon_get_standard_estimated_distance ( double *eta )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta);

   return fail;
}

int padcon_get_dobldobl_estimated_distance ( double *eta )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta);

   return fail;
}

int padcon_get_quaddobl_estimated_distance ( double *eta )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(887,&precision,b,eta);

   return fail;
}

int padcon_get_standard_hessian_step ( double *step )
{
   int fail;
   int precision = 0;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step);

   return fail;
}

int padcon_get_dobldobl_hessian_step ( double *step )
{
   int fail;
   int precision = 1;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step);

   return fail;
}

int padcon_get_quaddobl_hessian_step ( double *step )
{
   int fail;
   int precision = 2;
   int *b;

   fail = _ada_use_c2phc4c(888,&precision,b,step);

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

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(869,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

   *cre = nbr[0];
   *cim = nbr[1];

   return fail;
}

int padcon_get_quaddobl_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim )
{
   int fail;
   int inpars[3];
   double nbr[2];

   inpars[0] = 2;
   inpars[1] = 0;
   inpars[2] = leadidx;
   inpars[3] = idx;

   fail = _ada_use_c2phc4c(870,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr);

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

   fail = _ada_use_c2phc4c(871,inpars,&verbose,nbr);

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
