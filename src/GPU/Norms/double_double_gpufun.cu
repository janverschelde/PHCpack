// The file double_double_gpufun.cu defines the code for the functions
// specified in double_double_gpufun.h

#include <cmath>
#include "double_double_gpufun.h"

/************************** additions ********************************/

__device__ __forceinline__ double ddg_quick_two_sum
 ( double a, double b, double *err )
{
   double s = a + b;
   *err = b - (s - a);
   return s;
}

__device__ __forceinline__ double ddg_two_sum
 ( double a, double b, double *err )
{
   double s = a + b;
   double bb = s - a;
   *err = (a - (s - bb)) + (b - bb);
   return s;
}

__device__ __forceinline__ void ddg_add
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo )
{
   double s1, s2, t1, t2;

   s1 = ddg_two_sum(a_hi,b_hi,&s2);
   t1 = ddg_two_sum(a_lo,b_lo,&t2);
   s2 += t1;
   s1 = ddg_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   *c_hi = ddg_quick_two_sum(s1,s2,c_lo);
}

__device__ __forceinline__ double ddg_quick_two_diff
 ( double a, double b, double *err )
{
   double s = a - b;
   *err = (a - s) - b;
   return s;
}

__device__ __forceinline__ double ddg_two_diff
 ( double a, double b, double *err )
{
   double s = a - b;
   double bb = s - a;
   *err = (a - (s - bb)) - (b + bb);
   return s;
}

__device__ __forceinline__ void ddg_minus ( double *a_hi, double *a_lo )
{
   *a_hi = -(*a_hi);
   *a_lo = -(*a_lo);
}

__device__ __forceinline__ void ddg_sub
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo )
{
   double s1, s2, t1, t2;

   s1 = ddg_two_diff(a_hi,b_hi,&s2);
   t1 = ddg_two_diff(a_lo,b_lo,&t2);
   s2 += t1;
   s1 = ddg_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   *c_hi = ddg_quick_two_sum(s1,s2,c_lo);
}

__device__ __forceinline__ void ddg_sub_dd_d
 ( double a_hi, double a_lo, double b, double *c_hi, double *c_lo )
{
   double s1, s2;

   s1 = ddg_two_diff(a_hi,b,&s2);
   s2 += a_lo;
   *c_hi = ddg_quick_two_sum(s1,s2,c_lo);
}

/********** incrementers, decrementers, and multipliers ****************/

__device__ __forceinline__ void ddg_inc
 ( double *a_hi, double *a_lo, double b_hi, double b_lo )
{
   double s1, s2, t1, t2;

   s1 = ddg_two_sum(*a_hi,b_hi,&s2);
   t1 = ddg_two_sum(*a_lo,b_lo,&t2);
   s2 += t1;
   s1 = ddg_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   *a_hi = ddg_quick_two_sum(s1,s2,a_lo);
}

__device__ __forceinline__ void ddg_inc_d
 ( double *a_hi, double *a_lo, double b )
{
   double s1, s2;

   s1 = ddg_two_sum(*a_hi,b,&s2);
   s2 += *a_lo;
   *a_hi = ddg_quick_two_sum(s1,s2,a_lo);
}

__device__ __forceinline__ void ddg_dec
 ( double *a_hi, double *a_lo, double b_hi, double b_lo )
{
   double s1, s2, t1, t2;

   s1 = ddg_two_diff(*a_hi,b_hi,&s2);
   t1 = ddg_two_diff(*a_lo,b_lo,&t2);
   s2 += t1;
   s1 = ddg_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   *a_hi = ddg_quick_two_sum(s1,s2,a_lo);
}

__device__ __forceinline__ void ddg_dec_d
 ( double *a_hi, double *a_lo, double b )
{
   double s1, s2;

   s1 = ddg_two_diff(*a_hi,b,&s2);
   s2 += *a_lo;
   *a_hi = ddg_quick_two_sum(s1,s2,a_lo);
}

__device__ __forceinline__ void ddg_mlt
 ( double *a_hi, double *a_lo, double b_hi, double b_lo )
{
   double p1, p2;

   p1 = ddg_two_prod(*a_hi,b_hi,&p2);
   p2 += b_lo * (*a_lo);
   p2 += b_hi * (*a_hi);
   *a_hi = ddg_quick_two_sum(p1,p2,a_lo);
}

__device__ __forceinline__ void ddg_mlt_d
 ( double *a_hi, double *a_lo, double b )
{
   double p1, p2;

   p1 = ddg_two_prod(*a_hi,b,&p2);
   p2 += *a_lo * b;
   *a_hi = ddg_quick_two_sum(p1,p2,a_lo);
}

/************************ multiplications ********************************/

__device__ __forceinline__ void ddg_split ( double a, double *hi, double *lo )
{
   const double QD_SPLITTER = 134217729.0;            /* 2^27 + 1 */
   const double QD_SPLIT_THRESH = 6.69692879491417e+299; /* 2^996 */

   double temp;

   if (a > QD_SPLIT_THRESH || a < -QD_SPLIT_THRESH)
   {
      a *= 3.7252902984619140625e-09;  /* 2^-28 */
      temp = QD_SPLITTER * a;
      *hi = temp - (temp - a);
      *lo = a - *hi;
      *hi *= 268435456.0;  /* 2^28 */
      *lo *= 268435456.0;  /* 2^28 */
   }
   else
   {
      temp = QD_SPLITTER * a;
      *hi = temp - (temp - a);
      *lo = a - *hi;
   }
}

__device__ __forceinline__ double ddg_two_prod
 ( double a, double b, double *err )
{
   double a_hi,a_lo,b_hi,b_lo;
   double p = a*b;

   ddg_split(a,&a_hi,&a_lo);
   ddg_split(b,&b_hi,&b_lo);
   *err = ((a_hi*b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;

   return p;
}

__device__ __forceinline__ double ddg_two_sqr ( double a, double *err )
{
   double hi,lo;
   double q = a*a;

   ddg_split(a,&hi,&lo);
   *err = ((hi*hi - q) + 2.0*hi*lo) + lo*lo;

   return q;
}

__device__  __forceinline__ void ddg_mul
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo )
{
   double p1, p2;

   p1 = ddg_two_prod(a_hi,b_hi,&p2);
   p2 += (a_hi * b_lo + a_lo * b_hi);
   *c_hi = ddg_quick_two_sum(p1,p2,c_lo);
}

__device__  __forceinline__ void ddg_sqr
 ( double a_hi, double a_lo, double *b_hi, double *b_lo )
{
   double p1, p2;

   p1 = ddg_two_sqr(a_hi,&p2);
   p2 += 2.0 * a_lo * a_hi;
   p2 += a_lo * a_lo;
   *b_hi = ddg_quick_two_sum(p1,p2,b_lo);
}

__device__  __forceinline__ void ddg_mul_d_dd
 ( double a, double b_hi, double b_lo, double *c_hi, double *c_lo )
{
   double p1, p2;

   p1 = ddg_two_prod(a,b_hi,&p2);
   p2 += (b_lo * a);
   *c_hi = ddg_quick_two_sum(p1,p2,c_lo);
}

/*************************** divisions ***************************/

__device__  __forceinline__ void ddg_div
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo )
{
   double q1, q2, q3;
   double acc_hi, acc_lo;

   q1 = a_hi/b_hi;                             /* approximate quotient */
   ddg_mul_d_dd(q1,b_hi,b_lo,&acc_hi,&acc_lo); /* acc = q1*b */
   ddg_sub(a_hi,a_lo,acc_hi,acc_lo,c_hi,c_lo); /* c = a - q1 * b; */
   q2 = *c_hi/b_hi;
   ddg_mul_d_dd(q2,b_hi,b_lo,&acc_hi,&acc_lo); /* acc = q2*b */ 
   ddg_dec(c_hi,c_lo,acc_hi,acc_lo);           /* c -= (q2 * b); */
   q3 = *c_hi/b_hi;
   *c_hi = ddg_quick_two_sum(q1,q2,c_lo);
   ddg_inc_d(c_hi,c_lo,q3);                    /* c = ddg_real(q1, q2) + q3; */
}

/*************************** sqrt and abs ***************************/

__device__  __forceinline__ void ddg_sqrt
 ( double a_hi, double a_lo, double *b_hi, double *b_lo )
{
  /* Use Karp's trick: if x is an approximation to sqrt(a), then
       sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision. */
  
   if((a_hi == 0.0) && (a_lo == 0.0))
   {
      *b_hi = 0.0; *b_lo = 0.0;
   }
   else if(a_hi < 0.0)
   {
      *b_hi = -1.0; *b_lo = 0.0;
   }
   else
   {
      const double x = 1.0/sqrt(a_hi);
      double ax = a_hi*x;

      double sqax_hi,sqax_lo,y_hi,y_lo;
      
      ddg_sqr(ax,0.0,&sqax_hi,&sqax_lo);

      ddg_sub(a_hi,a_lo,sqax_hi,sqax_lo,&y_hi,&y_lo);     

      *b_hi = ddg_two_sum(ax,y_hi*x*0.5,b_lo);
   }
}

__device__  __forceinline__ void ddg_abs
 ( double a_hi, double a_lo, double *b_hi, double *b_lo )
{
   if(a_hi < 0.0)
   {
      *b_hi = -a_hi;
      *b_lo = -a_lo;
   }
   else
   {
      *b_hi = a_hi;
      *b_lo = a_lo;
   }
}
