// The file quad_double_functions.cpp defines the code for the functions
// specified in quad_double_functions.h

#include <cmath>
#include <iostream>
#include <iomanip>
#include "double_double_functions.h"
#include "quad_double_functions.h"

/************************** renormalizations **************************/

void qdf_renorm4
 ( double f0, double f1, double f2, double f3, double f4,
   double *pr, double *r0, double *r1, double *r2, double *r3 )
{
   int ptr;

   if(f1 == 0.0)
   {
      *pr = f0;
      ptr = 0;
      *r0 = ddf_quick_two_sum(*pr,f2,pr);
   }
   else
   {
      *r0 = f0;
      *pr = f1;
      ptr = 1;
      *r1 = ddf_quick_two_sum(*pr,f2,pr);
   }
   if(*pr == 0.0)
   {
      if(ptr == 0)
         *pr = *r0;
      else
         *pr = *r1;
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = ddf_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f3,pr);
   }

   if(*pr == 0.0)
   {
      if(ptr == 0)
      {
         *pr = *r0;
      }
      else if(ptr == 1)
      {
         *pr = *r1;
      }
      else if(ptr == 2)
      {
         *pr = *r2;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f4,pr);
   }

   if(*pr == 0.0)
   {
      if(ptr == 0)
      {
         *pr = *r0;
      }
      else if(ptr == 1)
      {
         *pr = *r1;
      }
      else if(ptr == 2)
      {
         *pr = *r2;
      }
      else if(ptr == 3)
      {
         *pr = *r3;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if((ptr < 4) && (*pr != 0.0))
   {
      if(ptr == 0)
      {
         *r0 = *pr;
      }
      else if(ptr == 1)
      {
         *r1 = *pr;
      }
      else if(ptr == 2)
      {
         *r2 = *pr;
      }
      else if(ptr == 3)
      {
         *r3 = *pr;
      }
      ptr = ptr + 1;
   }
   if(ptr < 1)
   {
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0; *r0 = 0.0;
   }
   else if(ptr < 2)
   {
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0;
   }
   else if(ptr < 3)
   {
      *r3 = 0.0; *r2 = 0.0;
   }
   else if(ptr < 4)
   {
      *r3 = 0.0;
   }
}

void qdf_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4,
   double *r0, double *r1, double *r2, double *r3 )
{
   double f0,f1,f2,f3,f4,pr;

   pr = ddf_quick_two_sum(x3,x4,&f4);
   pr = ddf_quick_two_sum(x2,pr,&f3);
   pr = ddf_quick_two_sum(x1,pr,&f2);
   f0 = ddf_quick_two_sum(x0,pr,&f1);

   qdf_renorm4(f0,f1,f2,f3,f4,&pr,r0,r1,r2,r3);
}

void qdf_renorm_add1
 ( double x0, double x1, double x2, double x3, double y,
   double *r0, double *r1, double *r2, double *r3 )
{
   double f0,f1,f2,f3,f4,pr;

   pr = ddf_two_sum(x3,y,&f4);
   pr = ddf_two_sum(x2,pr,&f3);
   pr = ddf_two_sum(x1,pr,&f2);
   f0 = ddf_two_sum(x0,pr,&f1);

   qdf_renorm4(f0,f1,f2,f3,f4,&pr,r0,r1,r2,r3);
}

/************************ copy and abs *******************************/

void qdf_copy
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double *b_hihi, double *b_lohi, double *b_hilo, double *b_lolo )
{
   *b_hihi = a_hihi;
   *b_lohi = a_lohi;
   *b_hilo = a_hilo;
   *b_lolo = a_lolo;
}

void qdf_abs
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double *b_hihi, double *b_lohi, double *b_hilo, double *b_lolo )
{
   if(a_hihi < 0.0)
   {
      *b_hihi = -a_hihi;
      *b_lohi = -a_lohi;
      *b_hilo = -a_hilo;
      *b_lolo = -a_lolo;
   }
   else
   {
      *b_hihi = a_hihi;
      *b_lohi = a_lohi;
      *b_hilo = a_hilo;
      *b_lolo = a_lolo;
   }
}

/****************** additions and subtractions ************************/

void qdf_add
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   // ALGORITHM : baileyAdd_fast<4,4,4> generated by CAMPARY.

   double f0,f1,f2,f3,f4,e;

   f4 = 0.0;
   f3 = ddf_two_sum(a_lolo,b_lolo,&e);
   f4 += e;
   f2 = ddf_two_sum(a_hilo,b_hilo,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_sum(a_lohi,b_lohi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f0 = ddf_two_sum(a_hihi,b_hihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;

   qdf_fast_renorm(f0,f1,f2,f3,f4,c_hihi,c_lohi,c_hilo,c_lolo);
}

void qdf_inc
 ( double *a_hihi, double *a_lohi, double *a_hilo, double *a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo )
{
   double f0,f1,f2,f3,f4,e;

   f4 = 0.0;
   f3 = ddf_two_sum(*a_lolo,b_lolo,&e);
   f4 += e;
   f2 = ddf_two_sum(*a_hilo,b_hilo,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_sum(*a_lohi,b_lohi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f0 = ddf_two_sum(*a_hihi,b_hihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;

   qdf_fast_renorm(f0,f1,f2,f3,f4,a_hihi,a_lohi,a_hilo,a_lolo);
}

void qdf_inc_d
 ( double *a_hihi, double *a_lohi, double *a_hilo, double *a_lolo,
   double b )
{
   qdf_renorm_add1(*a_hihi,*a_lohi,*a_hilo,*a_lolo,b,
                   a_hihi,a_lohi,a_hilo,a_lolo);
}

void qdf_dec
 ( double *a_hihi, double *a_lohi, double *a_hilo, double *a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo )
{
   double mbhihi = -b_hihi;
   double mblohi = -b_lohi;
   double mbhilo = -b_hilo;
   double mblolo = -b_lolo;

   qdf_inc(a_hihi,a_lohi,a_hilo,a_lolo,
           mbhihi,mblohi,mbhilo,mblolo);
}

void qdf_minus
 ( double *a_hihi, double *a_lohi, double *a_hilo, double *a_lolo )
{
   *a_hihi = -(*a_hihi);
   *a_lohi = -(*a_lohi);
   *a_hilo = -(*a_hilo);
   *a_lolo = -(*a_lolo);
}

void qdf_sub
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   qdf_copy(b_hihi,b_lohi,b_hilo,b_lolo,c_hihi,c_lohi,c_hilo,c_lolo);
   qdf_minus(c_hihi,c_lohi,c_hilo,c_lolo);
   qdf_inc(c_hihi,c_lohi,c_hilo,c_lolo,a_hihi,a_lohi,a_hilo,a_lolo);
}

/***************** multiplications and division ********************/

void qdf_mul_pwr2
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo, double b,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   *c_hihi = a_hihi*b;
   *c_lohi = a_lohi*b;
   *c_hilo = a_hilo*b;
   *c_lolo = a_lolo*b;
}

void qdf_mul
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   // ALGORITHM : baileyMul_fast<4,4,4> generated by CAMPARY.

   double f0,f1,f2,f3,f4,p,e;

   f4 =  a_lohi*b_lolo;
   f4 += a_hilo*b_hilo;
   f4 += a_lolo*b_lohi;
   f3 = ddf_two_prod(a_hihi,b_lolo,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,b_hilo,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   p = ddf_two_prod(a_hilo,b_lohi,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   p = ddf_two_prod(a_lolo,b_hihi,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   f2 = ddf_two_prod(a_hihi,b_hilo,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,b_lohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_hilo,b_hihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_prod(a_hihi,b_lohi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,b_hihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f0 = ddf_two_prod(a_hihi,b_hihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;

   qdf_fast_renorm(f0,f1,f2,f3,f4,c_hihi,c_lohi,c_hilo,c_lolo);
}

void qdf_sqr
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   double f0,f1,f2,f3,f4,p,e;

   f4 =  a_lohi*a_lolo;
   f4 += a_hilo*a_hilo;
   f4 += a_lolo*a_lohi;
   f3 = ddf_two_prod(a_hihi,a_lolo,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,a_hilo,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   p = ddf_two_prod(a_hilo,a_lohi,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   p = ddf_two_prod(a_lolo,a_hihi,&e);
   f4 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 += e;
   f2 = ddf_two_prod(a_hihi,a_hilo,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,a_lohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_hilo,a_hihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_prod(a_hihi,a_lohi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   p = ddf_two_prod(a_lohi,a_hihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f0 = ddf_two_prod(a_hihi,a_hihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;

   qdf_fast_renorm(f0,f1,f2,f3,f4,c_hihi,c_lohi,c_hilo,c_lolo);
}

void qdf_mul_qd_d
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double b,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   // ALGORITHM : baileyMul_fast<4,1,4>

   double f0,f1,f2,f3,f4,e;

   f4 = 0.0;
   f3 = ddf_two_prod(a_lolo,b,&e);
   f4 += e;
   f2 = ddf_two_prod(a_hilo,b,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f1 = ddf_two_prod(a_lohi,b,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;
   f0 = ddf_two_prod(a_hihi,b,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 += e;

   qdf_fast_renorm(f0,f1,f2,f3,f4,c_hihi,c_lohi,c_hilo,c_lolo);
}

void qdf_mlt_d
 ( double *a_hihi, double *a_lohi, double *a_hilo, double *a_lolo,
   double b )
{
   double c_hihi,c_lohi,c_hilo,c_lolo;

   qdf_mul_qd_d(*a_hihi,*a_lohi,*a_hilo,*a_lolo,b,
                &c_hihi,&c_lohi,&c_hilo,&c_lolo);

   *a_hihi = c_hihi;
   *a_lohi = c_lohi;
   *a_hilo = c_hilo;
   *a_lolo = c_lolo;
}

void qdf_div
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double b_hihi, double b_lohi, double b_hilo, double b_lolo,
   double *c_hihi, double *c_lohi, double *c_hilo, double *c_lolo )
{
   double acc_hihi,acc_lohi,acc_hilo,acc_lolo;
   double q0,q1,q2,q3,q4;

   q0 = a_hihi/b_hihi;
   qdf_mul_qd_d(b_hihi,b_lohi,b_hilo,b_lolo,q0,
                &acc_hihi,&acc_lohi,&acc_hilo,&acc_lolo);
   qdf_sub(a_hihi,a_lohi,a_hilo,a_lolo,
           acc_hihi,acc_lohi,acc_hilo,acc_lolo,c_hihi,c_lohi,c_hilo,c_lolo);

   q1 = *c_hihi/b_hihi;
   qdf_mul_qd_d(b_hihi,b_lohi,b_hilo,b_lolo,q1,
                &acc_hihi,&acc_lohi,&acc_hilo,&acc_lolo);
   qdf_sub(*c_hihi,*c_lohi,*c_hilo,*c_lolo,
           acc_hihi,acc_lohi,acc_hilo,acc_lolo,c_hihi,c_lohi,c_hilo,c_lolo);

   q2 = *c_hihi/b_hihi;
   qdf_mul_qd_d(b_hihi,b_lohi,b_hilo,b_lolo,q2,
                &acc_hihi,&acc_lohi,&acc_hilo,&acc_lolo);
   qdf_sub(*c_hihi,*c_lohi,*c_hilo,*c_lolo,
           acc_hihi,acc_lohi,acc_hilo,acc_lolo,c_hihi,c_lohi,c_hilo,c_lolo);

   q3 = *c_hihi/b_hihi;
   qdf_mul_qd_d(b_hihi,b_lohi,b_hilo,b_lolo,q3,
                &acc_hihi,&acc_lohi,&acc_hilo,&acc_lolo);
   qdf_sub(*c_hihi,*c_lohi,*c_hilo,*c_lolo,
           acc_hihi,acc_lohi,acc_hilo,acc_lolo,c_hihi,c_lohi,c_hilo,c_lolo);

   q4 = *c_hihi/b_hihi;

   qdf_fast_renorm(q0,q1,q2,q3,q4,c_hihi,c_lohi,c_hilo,c_lolo);
}

/***************************** square root *****************************/

void qdf_sqrt
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo,
   double *b_hihi, double *b_lohi, double *b_hilo, double *b_lolo )
{
   double z_hihi,z_lohi,z_hilo,z_lolo;

   ddf_sqrt(a_hihi,a_lohi,b_hihi,b_lohi);
   qdf_sqr(*b_hihi,*b_lohi,0.0,0.0,&z_hihi,&z_lohi,&z_hilo,&z_lolo);
   qdf_inc(&z_hihi,&z_lohi,&z_hilo,&z_lolo,a_hihi,a_lohi,a_hilo,a_lolo);
   qdf_div(z_hihi,z_lohi,z_hilo,z_lolo,*b_hihi,*b_lohi,0.0,0.0,
           &z_hihi,&z_lohi,&z_hilo,&z_lolo);
   qdf_mul_pwr2(z_hihi,z_lohi,z_hilo,z_lolo,0.5,
                b_hihi,b_lohi,b_hilo,b_lolo);
}

/*************************** basic output ***************************/

void qdf_write_doubles
 ( double a_hihi, double a_lohi, double a_hilo, double a_lolo )
{
   std::cout << std::scientific << std::setprecision(16);
   std::cout << "  hihi : " << a_hihi;
   std::cout << "  lohi : " << a_lohi << std::endl;
   std::cout << "  hilo : " << a_hilo;
   std::cout << "  lolo : " << a_lolo << std::endl;
}
