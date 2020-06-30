/* simple test on making a condition report in double precision */

#include <stdio.h>
#include <string.h>
#include "solcon.h"
#include "phcpack.h"

int main ( int argc, char *argv[] )
{
   int fail,len;
   char tofile,verbose,nls,sm;
   char name[80];
   int cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus,vrb;
   int t_err[16];
   int t_rco[16];
   int t_res[16];

   adainit();

   fail = read_standard_start_system();
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();

   fail = solcon_number_of_standard_solutions(&len);

   printf("\nRead %d solutions from file.\n",len);

   printf("Output to file ? (y/n) ");
   scanf("%c",&tofile);
   scanf("%c",&nls);     // swallow newline symbol

   printf("Verbose ? (y/n) ");
   scanf("%c",&verbose);
   scanf("%c",&nls);     // swallow newline symbol

   if(verbose == 'y')
      vrb = 1;
   else
      vrb = 0;

   if (tofile == 'n')
      fail = standard_condition_report
               (4,1.0e-8,1.0e-8,1.0e-8,0,name,
                &cntfail,&cntreal,&cntcmplx,&cntregu,&cntsing,&cntclus,
                t_err,t_rco,t_res,vrb);
   else 
   {
      printf("Give the name of the output file : ");
      scanf("%s",name);
      int nb = strlen(name);

      fail = standard_condition_report
               (4,1.0e-8,1.0e-8,1.0e-8,nb,name,
                &cntfail,&cntreal,&cntcmplx,&cntregu,&cntsing,&cntclus,
                t_err,t_rco,t_res,vrb);
   }
   printf("number of regular solutions   : %d\n", cntregu);
   printf("number of singular solutions  : %d\n", cntsing);
   printf("number of real solutions      : %d\n", cntreal);
   printf("number of complex solutions   : %d\n", cntcmplx);
   printf("number of clustered solutions : %d\n", cntclus);
   printf("number of failures            : %d\n", cntfail);

   printf("Frequency tables for forward errors, residual, ");
   printf("and condition numbers :\n");
   sm = 0;
   for(int i=0; i<16; i++)
   {
      printf(" %d", t_err[i]);
      sm = sm + t_err[i];
   }
   printf(" : %d\n", sm);
   sm = 0;
   for(int i=0; i<16; i++)
   {
      printf(" %d", t_res[i]);
      sm = sm + t_res[i];
   }
   printf(" : %d\n", sm);
   sm = 0;
   for(int i=0; i<16; i++)
   {
      printf(" %d", t_rco[i]);
      sm = sm + t_rco[i];
   }
   printf(" : %d\n", sm);

   printf("Small forward errors and residuals to the right.\n");
   printf("Well conditioned solutions are counted to the left.\n");

   adafinal();

   return 0;
}
