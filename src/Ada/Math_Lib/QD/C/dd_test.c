/* Test on the operations declared in the double_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "double_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   A test on the basic functions of double_double consists in
 *   generate two random doubles a, b and then calling the operations. */

int add_sub_of_pi_e();
/*
 * DESCRIPTION :
 *   Test on pi + e - pi. */

int div_sqr_of_pi ( void );
/*
 * DESCRIPTION :
 *   Divides pi^2 by pi. */

int log_exp_of_pi ( void );
/*
 * DESCRIPTION :
 *   Applies log after exp to the double double representation of pi. */

int my_sqrt ( void );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root. */

int main ( void )
{
   basic_test();
   add_sub_of_pi_e();
   div_sqr_of_pi();
   log_exp_of_pi();
   my_sqrt();
}

int basic_test ( void )
{
   double a,b,s,e,hi,lo;
   static const double _nan;

   printf("\ntesting some basic operations ...\n");
   srand(time(NULL));
   a = ((double) rand())/RAND_MAX;
   b = ((double) rand())/RAND_MAX;
   printf("  a = %21.14e\n  b = %21.14e\n",a,b);

   printf("nearest integer to a : %21.14e\n",dd_nint_d(a));
   printf("a truncated to integer : %21.14e\n",dd_aint(a));
   
   s = dd_two_sum(a,b,&e);
   printf("a+b = %21.14e  err = %21.14e\n",s,e);

   s = dd_two_diff(a,b,&e);
   printf("a-b = %21.14e  err = %21.14e\n",s,e);

   dd_split(a,&hi,&lo);
   printf("  a = %21.14e\n hi = %21.14e  lo = %21.14e\n",a,hi,lo);

   s = dd_two_prod(a,b,&e);
   printf("a*b = %21.14e  err = %21.14e\n",s,e);

   return 0;
}

int add_sub_of_pi_e ( void )
{
   const double pi_hi = 3.141592653589793116e+00; /* pi hi word */
   const double pi_lo = 1.224646799147353207e-16; /* pi lo word */
   double dd_pi[2];
   const double e_hi = 2.718281828459045091e+00; /* e hi word */
   const double e_lo = 1.445646891729250158e-16; /* e lo word */
   double dd_e[2];
   double pi_and_e[2];
   double pi_and_e_minus_pi[2];

   printf("\ntesting pi + e - pi ...\n");

   dd_pi[0] = pi_hi; dd_e[0] = e_hi;
   dd_pi[1] = pi_lo; dd_e[1] = e_lo;

   printf("      e hi = %21.14e         e lo = %21.14e\n",
          dd_e[0],dd_e[1]);

   dd_add(dd_pi,dd_e,pi_and_e);
   printf("   pi+e hi = %21.14e      pi+e lo = %21.14e\n",
          pi_and_e[0],pi_and_e[1]);

   dd_sub(pi_and_e,dd_pi,pi_and_e_minus_pi);
   printf("pi+e-pi hi = %21.14e   pi+e-pi lo = %21.14e\n",
          pi_and_e_minus_pi[0],pi_and_e_minus_pi[1]);

   return 0;
}

int div_sqr_of_pi ( void )
{
   const double pi_hi = 3.141592653589793116e+00; /* pi hi word */
   const double pi_lo = 1.224646799147353207e-16; /* pi lo word */
   double dd_pi[2];
   double sqr_of_pi[2];
   double div_of_sqr_of_pi[2];

   printf("\ntesting sqr(pi)/pi ...\n");

   dd_pi[0] = pi_hi;
   dd_pi[1] = pi_lo;
   printf("     pi hi = %21.14e        pi lo = %21.14e\n",
          dd_pi[0],dd_pi[1]);

   dd_sqr(dd_pi,sqr_of_pi);
   printf("   pi^2 hi = %21.14e      pi^2 lo = %21.14e\n",
          sqr_of_pi[0],sqr_of_pi[1]);

   dd_div(sqr_of_pi,dd_pi,div_of_sqr_of_pi);
   printf("pi^2/pi hi = %21.14e   pi^2/pi lo = %21.14e\n",
          div_of_sqr_of_pi[0],div_of_sqr_of_pi[1]);

   return 0;
}

int log_exp_of_pi ( void )
{
   const double pi_hi = 3.141592653589793116e+00; /* pi hi word */
   const double pi_lo = 1.224646799147353207e-16; /* pi lo word */
   double dd_pi[2];
   double exp_of_pi[2];
   double log_of_exp_of_pi[2];

   printf("\ntesting log(exp(pi)) = pi ...\n");

   dd_pi[0] = pi_hi;
   dd_pi[1] = pi_lo;
   printf("          pi hi = %21.14e            pi lo = %21.14e\n",
          dd_pi[0],dd_pi[1]);

   dd_exp(dd_pi,exp_of_pi);
   printf("     exp(pi) hi = %21.14e       exp(pi) lo = %21.14e\n",
          exp_of_pi[0],exp_of_pi[1]);

   dd_log(exp_of_pi,log_of_exp_of_pi);
   printf("log(exp(pi)) hi = %21.14e  log(exp(pi)) lo = %21.14e\n",
          log_of_exp_of_pi[0],log_of_exp_of_pi[1]);

   return 0;
}

int my_sqrt ( void )
{
   double n[2],x[2],y[2],z[2],e[2],a[2];
   const int max_steps = 9;
   char sqrt2[34] = "1.4142135623730950488016887242097\0";
   int i;

   n[0] = 2.0; n[1] = 0.0; dd_copy(n,x);

   puts("\nrunning Newton's method for sqrt(2) ...");
   dd_read(sqrt2,y);
   printf("y hi = %22.14e  y lo = %22.14e\n",y[0],y[1]);
   printf("step 0: "); dd_write(x,32); printf("\n");
   for(i=1; i <= max_steps; i++)
   {
      dd_mul(x,x,z);        /* z = x*x */
      dd_inc(z,n);          /* z += n */
      dd_div(z,x,z);        /* z /= x */
      dd_mlt_d(z,0.5);      /* z *= 0.5 */
      printf("step %d: ",i); dd_write(z,32);
      dd_copy(z,x);
      dd_sub(x,y,e);
      dd_abs(e,a);
      printf("  error : "); dd_write(a,3); printf("\n");
   }

   return 0;
}
