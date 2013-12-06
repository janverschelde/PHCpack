#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dcmplx.h"

void arithmetic_test ( dcmplx z1, dcmplx z2 );
/* tests addition, subtraction, multiplication, and division on the
   complex numbers z1 and z2 */

void random_test ( void );
/* generates random complex numbers to test the arithmetical operations */

void interactive_test ( void );
/* reads complex numbers given by user to test the operations on
   complex numbers */

int main(void)
{
   dcmplx one,imu,z1,z2;
   int ans;

   one = create1(1.0);

   printf("One as complex number : ");
   writeln_dcmplx(one);

   imu = create2(0.0,1.0);

   printf("The imaginary unit : ");
   writeln_dcmplx(imu);

   srand(time(NULL));

   for(;;)
   {
      printf("\nChoose one of the following :\n");
      printf("  0. exit this program; or\n");
      printf("  1. perform tests on user-given numbers; or\n");
      printf("  2. test operations on random numbers.\n");
      printf("Type 0, 1, or 2 to choose : ");
      scanf("%d", &ans);
      switch(ans)
      {
         case 0: return 0;
         case 1: interactive_test(); break;
         case 2: random_test(); break;
         default: printf("Invalid choice.  Please try again...\n");
      }
   }

   return 0;
}

void random_test ( void )
{
   dcmplx z1,z2;

   z1 = random_dcmplx1();

   printf("\nA random complex number z1 = ");
   writeln_dcmplx(z1);
   printf("  its modulus is "); printf("%.15le\n", modulus(z1));

   z2 = random_dcmplx1();

   printf("\nA random complex number z2 = ");
   writeln_dcmplx(z2);
   printf("  its modulus is "); printf("%.15le\n", modulus(z2));

   arithmetic_test(z1,z2);
}

void interactive_test ( void )
{
   dcmplx z1,z2;

   printf("Give complex number : ");
   read_dcmplx(&z1);

   printf("-> your number z1 : "); writeln_dcmplx(z1);
   printf("  its modulus : "); printf("%.15le\n", modulus(z1));

   printf("Give complex number : ");
   read_dcmplx(&z2);

   printf("-> your number z2 : "); writeln_dcmplx(z2);
   printf("its modulus : "); printf("%.15le\n", modulus(z2));
   
   arithmetic_test(z1,z2);
}

void arithmetic_test ( dcmplx z1, dcmplx z2 )
{
   dcmplx sum, dif, prod, quot, result;

   sum = add_dcmplx(z1,z2);
   dif = sub_dcmplx(z1,z2);
   prod = mul_dcmplx(z1,z2);
   quot = div_dcmplx(z1,z2);

   printf("\n  z1 + z2 = "); writeln_dcmplx(sum);
   printf("\n  z1 - z2 = "); writeln_dcmplx(dif);
   printf("\n  z1 * z2 = "); writeln_dcmplx(prod);
   printf("\n  z1 / z2 = "); writeln_dcmplx(quot);

   printf("\nRestoring z1 from sum and product with z2 : \n");
   result = sub_dcmplx(sum,z2);
   printf("\n  z1 + z2 - z2 = "); writeln_dcmplx(result);
   result = add_dcmplx(dif,z2);
   printf("\n  z1 - z2 + z2 = "); writeln_dcmplx(result);
   result = div_dcmplx(prod,z2);
   printf("\n  z1 * z2 / z2 = "); writeln_dcmplx(result);
   result = mul_dcmplx(quot,z2);
   printf("\n  z1 / z2 * z2 = "); writeln_dcmplx(result);
}
