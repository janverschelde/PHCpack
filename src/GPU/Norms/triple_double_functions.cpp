// The file triple_double_functions.cpp defines the code for the functions
// specified in triple_double_functions.h

#include <cmath>
#include "double_double_functions.h"
#include "triple_double_functions.h"

/************************* normalizations ***************************/

void tdf_fast_renorm
 ( double x0, double x1, double x2, double x3,
   double *r0, double *r1, double *r2 )
{
   double f0,f1,f2,f3,pr;
   int ptr;

   pr = ddf_quick_two_sum(x2,x3,&f3);
   pr = ddf_quick_two_sum(x1,pr,&f2);
   f0 = ddf_quick_two_sum(x0,pr,&f1);
   if(f1 == 0.0)
   {
      pr = f0;
      ptr = 0;
      *r0 = ddf_quick_two_sum(pr,f2,&pr);
   }
   else
   {
      *r0 = f0;
      pr = f1;
      ptr = 1;
      *r1 = ddf_quick_two_sum(pr,f2,&pr);
   }
   if(pr == 0.0)
   {
      if(ptr == 0)
         pr = *r0;
      else
         pr = *r1;
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
      *r0 = ddf_quick_two_sum(pr,f3,&pr);
   else if(ptr == 1)
      *r1 = ddf_quick_two_sum(pr,f3,&pr);
   else
      *r2 = ddf_quick_two_sum(pr,f3,&pr);

   if(pr == 0.0) 
   {
      if(ptr == 0)
         pr = *r0;
      else if(ptr == 1)
         pr = *r1;
      else
         pr = *r2;
    }
    else
    {
       ptr = ptr + 1;
    }
    if((ptr < 3) && (pr != 0.0))
    {
       if(ptr == 0)
          *r0 = pr;
       else if(ptr == 1)
          *r1 = pr;
       else
          *r2 = pr;
       ptr = ptr + 1;
    }
    if(ptr < 1)
    {
       *r0 = 0.0; *r1 = 0.0; *r2 = 0.0;
    }
    else if(ptr < 2)
    {
       *r1 = 0.0; *r2 = 0.0;
    }
    else if(ptr < 3)
    {
       *r2 = 0.0;
    }
}

/*************************** copy and abs ***************************/

void tdf_copy
 ( double a_hi, double a_mi, double a_lo,
   double *b_hi, double *b_mi, double *b_lo )
{
   *b_hi = a_hi;
   *b_mi = a_mi;
   *b_lo = a_lo;
}

void tdf_abs
 ( double a_hi, double a_mi, double a_lo,
   double *b_hi, double *b_mi, double *b_lo )
{
   if(a_hi < 0.0)
   {
      *b_hi = -a_hi;
      *b_mi = -a_mi;
      *b_lo = -a_lo;
   }
   else
   {
      *b_hi = a_hi;
      *b_mi = a_mi;
      *b_lo = a_lo;
   }
}

/************************** additions ********************************/

void tdf_add
 ( double a_hi, double a_mi, double a_lo,
   double b_hi, double b_mi, double b_lo,
   double *c_hi, double *c_mi, double *c_lo )
{
   double f0,f1,f2,f3,e;

   f2 = ddf_two_sum(a_lo,b_lo,&f3);
   f1 = ddf_two_sum(a_mi,b_mi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f0 = ddf_two_sum(a_hi,b_hi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;

   tdf_fast_renorm(f0,f1,f2,f3,c_hi,c_mi,c_lo);
}

void tdf_inc
 ( double *a_hi, double *a_mi, double *a_lo,
   double b_hi, double b_mi, double b_lo )
{
   double f0,f1,f2,f3,e;

   f2 = ddf_two_sum(*a_lo,b_lo,&f3);
   f1 = ddf_two_sum(*a_mi,b_mi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f0 = ddf_two_sum(*a_hi,b_hi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;

   tdf_fast_renorm(f0,f1,f2,f3,a_hi,a_mi,a_lo);
}

void tdf_minus ( double *a_hi, double *a_mi, double *a_lo )
{
   *a_hi = -(*a_hi);
   *a_mi = -(*a_mi);
   *a_lo = -(*a_lo);
}

void tdf_sub
 ( double a_hi, double a_mi, double a_lo,
   double b_hi, double b_mi, double b_lo,
   double *c_hi, double *c_mi, double *c_lo )
{
   tdf_copy(b_hi,b_mi,b_lo,c_hi,c_mi,c_lo);
   tdf_minus(c_hi,c_mi,c_lo);
   tdf_inc(c_hi,c_mi,c_lo,a_hi,a_mi,a_lo);
}

/***************** multiplications and division ********************/

void tdf_mul
 ( double a_hi, double a_mi, double a_lo,
   double b_hi, double b_mi, double b_lo,
   double *c_hi, double *c_mi, double *c_lo )
{
   double f0,f1,f2,f3,p,e;

   f3 = a_mi*b_lo + a_lo*b_mi;
   f2 = ddf_two_prod(a_hi,b_lo,&e);
   f3 += e;
   p = ddf_two_prod(a_mi,b_mi,&e);
   f3 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 += e;
   p = ddf_two_prod(a_lo,b_hi,&e);
   f3 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 += e;
   f1 = ddf_two_prod(a_hi,b_mi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   p = ddf_two_prod(a_mi,b_hi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f0 = ddf_two_prod(a_hi,b_hi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;

   tdf_fast_renorm(f0,f1,f2,f3,c_hi,c_mi,c_lo);
}

void tdf_sqr
 ( double a_hi, double a_mi, double a_lo,
   double *c_hi, double *c_mi, double *c_lo )
{
   double f0,f1,f2,f3,p,e;

   f3 = 2.0*a_mi*a_lo;
   f2 = ddf_two_prod(a_hi,a_lo,&e);
   f3 += e;
   p = ddf_two_prod(a_mi,a_mi,&e);
   f3 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 += e;
   p = ddf_two_prod(a_lo,a_hi,&e);
   f3 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 += e;
   f1 = ddf_two_prod(a_hi,a_mi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   p = ddf_two_prod(a_mi,a_hi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f0 = ddf_two_prod(a_hi,a_hi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;

   tdf_fast_renorm(f0,f1,f2,f3,c_hi,c_mi,c_lo);
}

void tdf_mul_td_d
 ( double a_hi, double a_mi, double a_lo, double b,
   double *c_hi, double *c_mi, double *c_lo )
{
   double f0,f1,f2,f3,e;

   f3 = 0.0;
   f2 = ddf_two_prod(a_lo,b,&e);
   f3 += e;
   f1 = ddf_two_prod(a_mi,b,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;
   f0 = ddf_two_prod(a_hi,b,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 += e;

   tdf_fast_renorm(f0,f1,f2,f3,c_hi,c_mi,c_lo);
}

void tdf_div
 ( double a_hi, double a_mi, double a_lo,
   double b_hi, double b_mi, double b_lo,
   double *c_hi, double *c_mi, double *c_lo )
{
   double acc_hi,acc_mi,acc_lo;
   double q0,q1,q2,q3;

   q0 = a_hi/b_hi;
   tdf_mul_td_d(b_hi,b_mi,b_lo,q0,&acc_hi,&acc_mi,&acc_lo);
   tdf_sub(a_hi,a_mi,a_lo,acc_hi,acc_mi,acc_lo,c_hi,c_mi,c_lo);

   q1 = *c_hi/b_hi;
   tdf_mul_td_d(b_hi,b_mi,b_lo,q1,&acc_hi,&acc_mi,&acc_lo);
   tdf_sub(*c_hi,*c_mi,*c_lo,acc_hi,acc_mi,acc_lo,c_hi,c_mi,c_lo);

   q2 = *c_hi/b_hi;
   tdf_mul_td_d(b_hi,b_mi,b_lo,q2,&acc_hi,&acc_mi,&acc_lo);
   tdf_sub(*c_hi,*c_mi,*c_lo,acc_hi,acc_mi,acc_lo,c_hi,c_mi,c_lo);

   q3 = *c_hi/b_hi;

   tdf_fast_renorm(q0,q1,q2,q3,c_hi,c_mi,c_lo);
}
