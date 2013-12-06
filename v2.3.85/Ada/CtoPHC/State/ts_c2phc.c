/* simple test program to the operations in PHCpack */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_c_to_phcpack ( int n );
extern void adafinal();

int main(void)
{
   printf("\nTesting C to PHCpack...\n");

   adainit();
   {
      int fail,choice = 0;
      char skip_newline;
      do
      { 
         fail = _ada_c_to_phcpack(choice);
         printf("Give a number (0 for menu, -1 to exit) : ");
         scanf("%d",&choice);
         skip_newline = getchar();
      } while (choice >= 0);
   }
   adafinal();

   return 0;
}
