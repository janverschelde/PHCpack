/* simple test on running a Pade continuation */

#include <stdio.h>
#include "padcon.h"

void write_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Writes the current values of the homotopy continuation parameters. */

void tune_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Interactive loop to tune the homotopy continuation parameters. */

int main ( void )
{
   adainit();

   printf("Tuning the homotopy continuation parameters ...\n");

   padcon_set_default_parameters();
   tune_homotopy_continuation_parameters();

   adafinal();

   return 0;
}

void write_homotopy_continuation_parameters ( void )
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
   printf(" 6. multiplication factor of the series step   : ");
   fail = padcon_get_homotopy_continuation_parameter(6,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 7. multiplication factor of the pole radius   : ");
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

void tune_homotopy_continuation_parameters ( void )
{
   int choice,intpar,fail;
   double fltpar;
   double gamma[2];

   do
   {
      write_homotopy_continuation_parameters();
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
         printf("-> give a new multiplication factor for the series step : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 7)
      {
         printf("-> give a new multiplication factor for the pole radius : ");
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
