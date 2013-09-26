/* simple test program to the interactive operations of use_c2phc,
 * i.e.: this program is limited to those operations which require
 *       explicit dialogue with the user */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_c2phc ( int job, int *a, int *b, double *c );
extern void adafinal();

int main(void)
{
   printf("\nTesting using PHCpack from within C ...\n");

   adainit();
   {
      int fail,choice = 0;
      int *a,*b;
      double *c;
      char skip_newline;
      do
      { 
         fail = _ada_use_c2phc(choice,a,b,c);
         printf("Give a number (0 for menu, -1 to exit) : ");
         scanf("%d",&choice);
         skip_newline = getchar();
      } while (choice >= 0);
   }
   adafinal();

   return 0;
}
