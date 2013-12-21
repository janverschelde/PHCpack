/* simple test program to the interactive operations of use_c2fac,
 * i.e.: this program is limited to those operations which require
 *       explicit dialogue with the user */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_c2fac ( int job, int *a, int *b, double *c );
extern void adafinal();

int main(void)
{
   printf("\nTesting using factorization in PHCpack from within C ...\n");

   adainit();
   {
      int fail,choice = 0;
      int *a;
      int b[2];
      double *c;
      char skip_newline;
      do
      { 
         if(choice == 3)
         {
            double gamma[2];
            gamma[0] = 0.929016078002768;
            gamma[1] = 0.370039358463874;
            fail = _ada_use_c2fac(choice,a,b,gamma);
         }
         else if(choice == 5)
         {
            int n;
            printf("Give max #loops : "); scanf("%d",&n);
            printf("initialize with %d loops...\n",n);
            fail = _ada_use_c2fac(choice,&n,b,c);
         }
         else
         {
            fail = _ada_use_c2fac(choice,a,b,c);
            if(choice == 1) printf("fail = %d, n = %d, dim = %d, deg = %d\n",
                                   fail,*a,b[0],b[1]);
         }
         printf("Give a number (0 for menu, -1 to exit) : ");
         scanf("%d",&choice);
         skip_newline = getchar();
      } while (choice >= 0);
   }
   adafinal();

   return 0;
}
