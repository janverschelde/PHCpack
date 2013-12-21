/* simple test program in C on the Ada procedure use_syspool */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_syspool ( int task, int *a, int *b, double *c );
extern void adafinal();

int main(void)
{
   printf("\nTesting Systems Pool...\n");

   adainit();
   {
      int n,s,k,fail,*b;
      double *c;

      printf("\nGive the size of the pool : ");
      scanf("%d",&n);

      fail = _ada_use_syspool(0,&n,b,c);   /* initialize pool size*/
      fail = _ada_use_syspool(1,&s,b,c);   /* get the pool size */

      printf("The size of the pool after initialization : %d\n",s);
      printf("Give k (<= %d) : ",s); scanf("%d",&k); 
      fail = _ada_use_syspool(2,&k,b,c);   /* read k-th system */
      printf("System %d in the pool :\n",k); 
      fail = _ada_use_syspool(3,&k,b,c);   /* write k-th system */
   }
   adafinal();

   return 0;
}
