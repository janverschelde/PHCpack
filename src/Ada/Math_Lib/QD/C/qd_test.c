/* Test on the operations declared in the quad_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "quad_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   A test on the basic functions of quad_double. */

int add_sub_of_pi_e ( void );
/*
 * DESCRIPTION :
 *   Test on p + e - pi. */

int div_sqr_of_pi ( void );
/*
 * DESCRIPTION :
 *   Test on p^2/pi = pi. */

int log_exp_of_pi ( void );
/*
 * DESCRIPTION :
 *   Applies log after exp to the quad double representation of pi. */

int my_sqrt ( void );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root. */

int qd_write_pi ( void );
/*
 * DESCRIPTION :
 *   Writes the quad double precision value for Pi via a string. */

int main ( void )
{
   basic_test();
   add_sub_of_pi_e();
   div_sqr_of_pi();
   log_exp_of_pi();
   my_sqrt();
   qd_write_pi();
}

int basic_test ( void )
{
   double c0,c1,c2,c3,c4;

   printf("\ntesting some basic operations ...\n");
   srand(time(NULL));
   c0 = ((double) rand())/RAND_MAX;
   c1 = ((double) rand())/RAND_MAX;
   c2 = ((double) rand())/RAND_MAX;
   c3 = ((double) rand())/RAND_MAX;
   c4 = ((double) rand())/RAND_MAX;

   qd_quick_renorm(&c0,&c1,&c2,&c3,&c4);

   return 0;
}

int add_sub_of_pi_e ( void )
{
   const double qd_pi_0 =  3.141592653589793116e+00; /* qd_pi[0] */
   const double qd_pi_1 =  1.224646799147353207e-16; /* qd_pi[1] */
   const double qd_pi_2 = -2.994769809718339666e-33; /* qd_pi[2] */
   const double qd_pi_3 =  1.112454220863365282e-49; /* qd_pi[3] */
   double qd_pi[4];
   const double qd_e_0 =  2.718281828459045091e+00; /* qd_e[0] */
   const double qd_e_1 =  1.445646891729250158e-16; /* qd_e[1] */
   const double qd_e_2 = -2.127717108038176765e-33; /* qd_e[2] */
   const double qd_e_3 =  1.515630159841218954e-49; /* qd_e[3] */
   double qd_e[4];
   double pi_and_e[4];
   double pi_and_e_minus_pi[4];

   printf("\ntesting pi + e - pi ...\n");

   qd_pi[0] = qd_pi_0; qd_pi[1] = qd_pi_1;
   qd_pi[2] = qd_pi_2; qd_pi[3] = qd_pi_3;
   qd_e[0] = qd_e_0; qd_e[1] = qd_e_1;
   qd_e[2] = qd_e_2; qd_e[3] = qd_e_3;

   printf("pi[0] = %22.14e  e[0] = %22.14e\n", qd_pi[0],qd_e[0]);
   printf("pi[1] = %22.14e  e[1] = %22.14e\n", qd_pi[1],qd_e[1]);
   printf("pi[2] = %22.14e  e[2] = %22.14e\n", qd_pi[2],qd_e[2]);
   printf("pi[3] = %22.14e  e[3] = %22.14e\n", qd_pi[3],qd_e[3]);

   qd_add(qd_pi,qd_e,pi_and_e);
   qd_sub(pi_and_e,qd_pi,pi_and_e_minus_pi);

   printf("pi+e[0] = %22.14e  pi+e-pi[0] = %22.14e\n",
          pi_and_e[0],pi_and_e_minus_pi[0]);
   printf("pi+e[1] = %22.14e  pi+e-pi[1] = %22.14e\n",
          pi_and_e[1],pi_and_e_minus_pi[1]);
   printf("pi+e[2] = %22.14e  pi+e-pi[2] = %22.14e\n",
          pi_and_e[2],pi_and_e_minus_pi[2]);
   printf("pi+e[3] = %22.14e  pi+e-pi[3] = %22.14e\n",
          pi_and_e[3],pi_and_e_minus_pi[3]);

   return 0;
}

int div_sqr_of_pi ( void )
{
   const double qd_pi_0 =  3.141592653589793116e+00; /* qd_pi[0] */
   const double qd_pi_1 =  1.224646799147353207e-16; /* qd_pi[1] */
   const double qd_pi_2 = -2.994769809718339666e-33; /* qd_pi[2] */
   const double qd_pi_3 =  1.112454220863365282e-49; /* qd_pi[3] */
   double qd_pi[4];
   double sqr_pi[4];
   double pi_times_pi[4];
   double div_of_sqr_of_pi[4];

   qd_pi[0] = qd_pi_0; qd_pi[1] = qd_pi_1;
   qd_pi[2] = qd_pi_2; qd_pi[3] = qd_pi_3;

   printf("\ntesting pi^2/pi ...\n");

   qd_sqr(qd_pi,sqr_pi);
   qd_mul(qd_pi,qd_pi,pi_times_pi);
   qd_div(sqr_pi,qd_pi,div_of_sqr_of_pi);

   printf("pi[0] = %22.14e  pi^2[0] = %22.14e\n",qd_pi[0],sqr_pi[0]);
   printf("pi[1] = %22.14e  pi^2[1] = %22.14e\n",qd_pi[1],sqr_pi[1]);
   printf("pi[2] = %22.14e  pi^2[2] = %22.14e\n",qd_pi[2],sqr_pi[2]);
   printf("pi[3] = %22.14e  pi^2[3] = %22.14e\n",qd_pi[3],sqr_pi[3]);

   printf("pi*pi[0] = %22.14e  pi^2/pi[0] = %22.14e\n",
          pi_times_pi[0],div_of_sqr_of_pi[0]);
   printf("pi*pi[1] = %22.14e  pi^2/pi[1] = %22.14e\n",
          pi_times_pi[1],div_of_sqr_of_pi[1]);
   printf("pi*pi[2] = %22.14e  pi^2/pi[2] = %22.14e\n",
          pi_times_pi[2],div_of_sqr_of_pi[2]);
   printf("pi*pi[3] = %22.14e  pi^2/pi[3] = %22.14e\n",
          pi_times_pi[3],div_of_sqr_of_pi[3]);

   return 0;
}

int log_exp_of_pi ( void )
{
   const double qd_pi_0 =  3.141592653589793116e+00; /* qd_pi[0] */
   const double qd_pi_1 =  1.224646799147353207e-16; /* qd_pi[1] */
   const double qd_pi_2 = -2.994769809718339666e-33; /* qd_pi[2] */
   const double qd_pi_3 =  1.112454220863365282e-49; /* qd_pi[3] */
   double qd_pi[4];
   double exp_of_pi[4];
   double log_of_exp_of_pi[4];

   qd_pi[0] = qd_pi_0; qd_pi[1] = qd_pi_1;
   qd_pi[2] = qd_pi_2; qd_pi[3] = qd_pi_3;

   printf("\ntesting log(exp(pi)) = pi ...\n");

   qd_exp(qd_pi,exp_of_pi);
   printf("pi[0] = %22.14e  exp(pi)[0] = %22.14e\n", qd_pi[0],exp_of_pi[0]);
   printf("pi[1] = %22.14e  exp(pi)[1] = %22.14e\n", qd_pi[1],exp_of_pi[1]);
   printf("pi[2] = %22.14e  exp(pi)[2] = %22.14e\n", qd_pi[2],exp_of_pi[2]);
   printf("pi[3] = %22.14e  exp(pi)[3] = %22.14e\n", qd_pi[3],exp_of_pi[3]);

   qd_log(exp_of_pi,log_of_exp_of_pi);
   printf("pi[0] = %22.14e  log(exp(pi))[0] = %22.14e\n",
          qd_pi[0],log_of_exp_of_pi[0]);
   printf("pi[1] = %22.14e  log(exp(pi))[1] = %22.14e\n",
          qd_pi[1],log_of_exp_of_pi[1]);
   printf("pi[2] = %22.14e  log(exp(pi))[2] = %22.14e\n",
          qd_pi[2],log_of_exp_of_pi[2]);
   printf("pi[3] = %22.14e  log(exp(pi))[3] = %22.14e\n",
          qd_pi[3],log_of_exp_of_pi[3]);

   return 0;
}

int my_sqrt ( void )
{
   double n[4],x[4],y[4],z[4],e[4],a[4];
   const int max_steps = 9;
   char sqrt2[66]
      = "1.414213562373095048801688724209698078569671875376948073176679738\0";
   int i;

   qd_real(2.0,n);
   qd_copy(n,x);

   puts("\nrunning Newton's method for sqrt(2) ...");
   printf("applying qd_read to the string \n%s\n",sqrt2);
   qd_read(sqrt2,y);
   for(i=0; i<4; i++) printf("y[%d] = %22.14e\n",i,y[i]);

   printf(" sqrt2: "); qd_write(y,64); printf("\n");
   printf("step 0: "); qd_write(x,64); printf("\n");
   for(i=1; i <= max_steps; i++)
   {
      qd_mul(x,x,z);        /* z = x*x */
      qd_add(z,n,z);        /* z += n */
      qd_div(z,x,z);        /* z /= x */
      qd_mul_d(z,0.5,z);    /* z *= 0.5 */
      printf("step %d: ",i); qd_write(z,64);
      qd_copy(z,x);
      qd_sub(x,y,e);
      qd_abs(e,a);
      printf("  error : "); qd_write(a,3); printf("\n");
   }
   return 0;
}

int qd_write_pi ( void )
{
   char s[74];
   int s_end;
   const double qd_pi_0 =  3.141592653589793116e+00; /* qd_pi[0] */
   const double qd_pi_1 =  1.224646799147353207e-16; /* qd_pi[1] */
   const double qd_pi_2 = -2.994769809718339666e-33; /* qd_pi[2] */
   const double qd_pi_3 =  1.112454220863365282e-49; /* qd_pi[3] */
   double qd_pi[4];
   qd_pi[0] = qd_pi_0;
   qd_pi[1] = qd_pi_1;
   qd_pi[2] = qd_pi_2;
   qd_pi[3] = qd_pi_3;

   printf("writing Pi with quad double precision ...\n");
   qd_to_string(qd_pi,s,64,0,0,0,1,' ',&s_end);
   printf("Pi = %s\n",s);

   return 0;
}
