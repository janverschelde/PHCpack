/* testing the determination of the continuation parameters */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"

void write_continuation_parameters ( double *v );
/* writes the current values for the continuation parameters to screen */

int main ( int argc, char *argv[] )
{
   int fail,i;
   double v[34],x;

   adainit();

   printf("\nTuning the Continuation parameters in PHCpack...\n");
   /* fail = tune_continuation_parameters(); */
   printf("\n");
   do
   {
      fail = retrieve_continuation_parameters(v);
      write_continuation_parameters(v);
      printf("Give index of parameter to change (0 to quit) : ");
      scanf("%d",&i);
      if (i > 0)
      {
         printf("Give value for parameter %d : ",i); scanf("%lf",&x);
         v[i-1] = x;
         set_continuation_parameters(v);
         {
            char ans = 'y';
            do
            {
               fail = get_value_of_continuation_parameter(i,&x);
               printf("-> the value for parameter %d : %.3e.\n",i,x);
               printf("-> change value of index ? (y/n) ");
               scanf("%d",&ans); /* skip previous newline symbol */
               scanf("%c",&ans);
               if(ans != 'y') break;
               scanf("%c",&ans); /* skip the new line */
               printf("-> give value for index %d : ",i);
               scanf("%lf",&x);
               fail = set_value_of_continuation_parameter(i,&x);
            }
            while (ans != 0);
         }
      }
   } while ( i > 0);

   adafinal();

   return 0;
}

const char *predictor_banner ( int pred_type )
/* returns banner for corresponding predictor type */
{
   switch (pred_type)
   {
      case 0:  return "( x:Sec,t:Rea )";
      case 1:  return "( x:Sec,t:Com )"; 
      case 2:  return "( x:Sec,t:Geo )";
      case 3:  return "( x:Tan,t:Rea )";
      case 4:  return "( x:Tan,t:Com )";
      case 5:  return "( x:Tan,t:Geo )";
      case 6:  return "( x:Her,t:Rea )";
      default: return " no predictor  ";
   }
}

void write_continuation_parameters ( double *v )
{
   printf
   ("****************** CURRENT CONTINUATION PARAMETERS *****************\n");
   printf("GLOBAL MONITOR : \n");
   printf("  1. the condition of the homotopy           : %d\n",(int)v[0]);
   printf("  2. number of paths tracked simultaneously  : %d\n",(int)v[1]);
   printf("  3. maximum number of steps along a path    : %d\n",(int)v[2]); 
   printf("  4. distance from target to start end game  : %.3e\n",v[3]); 
   printf("  5. order of extrapolator in end game       : %d\n",(int)v[4]);
   printf("  6. maximum number of re-runs               : %d\n",(int)v[5]);
   printf("STEP CONTROL (PREDICTOR) : ");
   printf("                   along path : end game\n");
   printf("  7: 8. type %s:%s : ",
          predictor_banner((int)v[6]),predictor_banner((int)v[7]));
   printf("%d",(int)v[6]);             printf("         : ");
   printf("%d\n",(int)v[7]);
   printf("  9:10. minimum step size                    : ");
   printf("%.3e",v[8]);                        printf(" : ");
   printf("%.3e\n",v[9]);
   printf(" 11:12. maximum step size                    : ");
   printf("%.3e",v[10]);                       printf(" : ");
   printf("%.3e\n",v[11]);
   printf(" 13:14. reduction factor for step size       : ");
   printf("%.3e",v[12]);                       printf(" : ");
   printf("%.3e\n",v[13]);
   printf(" 15:16. expansion factor for step size       : ");
   printf("%.3e",v[15]);                       printf(" : ");
   printf("%.3e\n",v[16]);
   printf(" 17:18. expansion threshold                  : ");
   printf("%d",(int)v[16]);            printf("         : ");
   printf("%d\n",(int)v[17]);
   printf("PATH CLOSENESS (CORRECTOR) : ");
   printf("                 along path : end game\n");
   printf(" 19:20. maximum number of iterations         : ");
   printf("%d",(int)v[18]);            printf("         : ");
   printf("%d\n",(int)v[19]);
   printf(" 21:22. relative precision for residuals     : ");
   printf("%.3e",v[20]);                       printf(" : ");
   printf("%.3e\n",v[21]); 
   printf(" 23:24. absolute precision for residuals     : ");
   printf("%.3e",v[22]);                       printf(" : ");
   printf("%.3e\n",v[23]);
   printf(" 25:26. relative precision for corrections   : ");
   printf("%.3e",v[24]);                       printf(" : ");
   printf("%.3e\n",v[25]);
   printf(" 27:28. absolute precision for corrections   : ");
   printf("%.3e",v[26]);                       printf(" : ");
   printf("%.3e\n",v[27]); 
   printf("SOLUTION TOLERANCES : ");
   printf("                        along path : end game\n");
   printf(" 29:30. inverse condition of Jacobian        : ");
   printf("%.3e",v[28]);                       printf(" : ");
   printf("%.3e\n",v[29]);
   printf(" 31:32. clustering of solutions              : ");
   printf("%.3e",v[30]);                       printf(" : ");
   printf("%.3e\n",v[31]);
   printf(" 33:34. solution at infinity                 : ");
   printf("%.3e",v[32]);                       printf(" : ");
   printf("%.3e\n",v[33]);
   printf
   ("********************************************************************\n");
}
