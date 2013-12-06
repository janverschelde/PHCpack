/* simple test program in C on the Ada procedure use_solpool */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_solpool ( int task, int *a, int *b, double *c );
extern void adafinal();

int main(void)
{
   printf("\nTesting Solutions Pool...\n");

   adainit();
   {
      int n,s,fail,*b;
      double *c;

      printf("\nGive the size of the pool : ");
      scanf("%d",&n);

      fail = _ada_use_solpool(0,&n,b,c);   /* initialize pool size*/
      fail = _ada_use_solpool(1,&s,b,c);   /* get the pool size */

      printf("The size of the pool after initialization : %d\n",s);
   }
   adafinal();

   return 0;
}
