/* Test of "time.c" on computing sums and products 
   of natural numbers. Do "time ts_timer" when running 
   this program to see whether the UNIX timer gives the 
   same results as the printed user CPU times. 
   User can also choose a large number in milliseconds to 
   test if the printed time format is correct.  */

#include <stdio.h>
#include "timer.h"
#include <time.h>

int main()
{
  timer total, local;
  int acc, n, m, choice, i;
  long long mseconds;

  tstart(&total);

  printf("Please choose a test method:\n");
  printf("  1 Do time ts_timer to see if the UNIX timer gives");
  printf(" the same results as the printed times.\n");
  printf("  2 Give the number of seconds to check if the printed");
  printf(" time format is correct\n");
  printf("Your choice is: ");
  scanf("%d", &choice);
  if(choice==1)
  {
     printf("Test of timer on computation of sums and products.\n");
     printf("Give the number of sums     :"); scanf("%d", &n);
     printf("Give the number of products :"); scanf("%d", &m);
     tstart(&local);
     acc = 0;
     for(i=1; i<=n; i++)
        acc = acc+1;
     printf("1 + 1 + ... + 1 : %d\n", acc);
     printf("Calculating the sum\n");
     tstop(&local);
     printf("Sum ");
     tprint(local);
     
     tstart(&local);
     acc = 1;
     for(i=1; i<=m; i++)
        acc = acc*1;     
     printf("1 * 1 * ... * 1 : %d\n", acc);
     printf("Calculating the product\n");
     tstop(&local);
     printf("Product ");
     tprint(local);
     tstop(&total);
     tprint(total);  
   }
   else if(choice==2)
   {
     printf("Please give the number of milliseconds to test: ");
     scanf("%lld", &mseconds);
     total.user_time = mseconds*CLOCKS_PER_SEC/1000;
     total.time1 = 0;
     total.time2 = mseconds/1000;
     printf("Total elapsed time is %.3lf seconds\n", (double)mseconds/1000);
     printf("Calculated ");
     tprint(total);
   }
   else
   {
     printf("Please choose 1 or 2.\n");
   }
        
  return 0;
}
