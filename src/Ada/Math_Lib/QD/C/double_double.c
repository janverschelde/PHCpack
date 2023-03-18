/* file: double_double.c */

/* This file contains the corresponding C code for the functions
   with prototypes declared in the double_double.h file. */

/* basic functions */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "double_double.h"

/* Part I: basic functions from inline.h */

double dd_quick_two_sum ( double a, double b, double *err )
{
   double s = a + b;
   *err = b - (s - a);
   return s;
}

double dd_quick_two_diff ( double a, double b, double *err )
{
   double s = a - b;
   *err = (a - s) - b;
   return s;
}

double dd_two_sum ( double a, double b, double *err )
{
   double s = a + b;
   double bb = s - a;
   *err = (a - (s - bb)) + (b - bb);
   return s;
}

double dd_two_diff ( double a, double b, double *err )
{
   double s = a - b;
   double bb = s - a;
   *err = (a - (s - bb)) - (b + bb);
   return s;
}

void dd_split ( double a, double *hi, double *lo )
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

double dd_two_prod ( double a, double b, double *err )
{
   double a_hi,a_lo,b_hi,b_lo;
   double p = a*b;

   dd_split(a,&a_hi,&a_lo);
   dd_split(b,&b_hi,&b_lo);
   *err = ((a_hi*b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;

   return p;
}

double dd_two_sqr ( double a, double *err )
{
   double hi,lo;
   double q = a*a;

   dd_split(a,&hi,&lo);
   *err = ((hi*hi - q) + 2.0*hi*lo) + lo*lo;

   return q;
}

double dd_nint_d ( double d )
{
   if (d == floor(d)) return d;
   return floor(d + 0.5);
}

double dd_aint ( double d )
{
   return (d >= 0.0) ? floor(d) : ceil(d);
}

/* Part II: arithmetical operations from dd_inline.h */

/*************************** additions ***************************/

void dd_add ( const double *a, const double *b, double *c )
{
/* satisfies IEEE style error bound, due to K. Briggs and W. Kahan */

   double s1, s2, t1, t2;

   s1 = dd_two_sum(a[0],b[0],&s2);
   t1 = dd_two_sum(a[1],b[1],&t2);
   s2 += t1;
   s1 = dd_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_add_d_dd ( double a, const double *b, double *c )
{
   double s1, s2;

   s1 = dd_two_sum(b[0],a,&s2);
   s2 += b[1];
   c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_add_dd_d ( const double *a, double b, double *c )
{
   double s1, s2;

   s1 = dd_two_sum(a[0],b,&s2);
   s2 += a[1];
   c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_add_d_d ( double a, double b, double *c )
{
   c[0] = dd_two_sum(a,b,&c[1]);
}

void dd_inc ( double *a, const double *b )
{
/* we assume that QD_IEEE_ADD is defined, see dd_add comment */

  double s1, s2, t1, t2;

  s1 = dd_two_sum(a[0],b[0],&s2);
  t1 = dd_two_sum(a[1],b[1],&t2);
  s2 += t1;
  s1 = dd_quick_two_sum(s1,s2,&s2);
  s2 += t2;
  a[0] = dd_quick_two_sum(s1,s2,&a[1]);
}

void dd_inc_d ( double *a, double b )
{
   double s1, s2;

   s1 = dd_two_sum(a[0],b,&s2);
   s2 += a[1];
   a[0] = dd_quick_two_sum(s1,s2,&a[1]);
}

/************************* subtractions **************************/

void dd_sub ( const double *a, const double *b, double *c )
{
/* we assume that QD_IEEE_ADD is defined, see dd_add comment */

  double s1, s2, t1, t2;

  s1 = dd_two_diff(a[0],b[0],&s2);
  t1 = dd_two_diff(a[1],b[1],&t2);
  s2 += t1;
  s1 = dd_quick_two_sum(s1,s2,&s2);
  s2 += t2;
  c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_sub_d_dd ( double a, const double *b, double *c )
{
   double s1, s2;

   s1 = dd_two_diff(a,b[0],&s2);
   s2 -= b[1];
   c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_sub_dd_d ( const double *a, double b, double *c )
{
   double s1, s2;

   s1 = dd_two_diff(a[0],b,&s2);
   s2 += a[1];
   c[0] = dd_quick_two_sum(s1,s2,&c[1]);
}

void dd_sub_d_d ( double a, double b, double *c )
{
   c[0] = dd_two_diff(a,b,&c[1]);
}

void dd_dec ( double *a, const double *b )
{
/* we assume QD_IEEE_ADD is defined, see dd_add comment */

   double s1, s2, t1, t2;

   s1 = dd_two_diff(a[0],b[0],&s2);
   t1 = dd_two_diff(a[1],b[1],&t2);
   s2 += t1;
   s1 = dd_quick_two_sum(s1,s2,&s2);
   s2 += t2;
   a[0] = dd_quick_two_sum(s1,s2,&a[1]);
}

void dd_dec_d ( double *a, double b )
{
   double s1, s2;

   s1 = dd_two_diff(a[0],b,&s2);
   s2 += a[1];
   a[0] = dd_quick_two_sum(s1,s2,&a[1]);
}

/******************************* unary minus **************************/

void dd_minus ( double *a )
{
   a[0] = -a[0];
   a[1] = -a[1];
}

/*************************** multiplications ***************************/

void dd_mul ( const double *a, const double *b, double *c )
{
   double p1, p2;

   p1 = dd_two_prod(a[0],b[0],&p2);
   p2 += (a[0] * b[1] + a[1] * b[0]);
   c[0] = dd_quick_two_sum(p1,p2,&c[1]);
}

void dd_mul_d_dd ( double a, const double *b, double *c )
{
   double p1, p2;

   p1 = dd_two_prod(a,b[0],&p2);
   p2 += (b[1] * a);
   c[0] = dd_quick_two_sum(p1,p2,&c[1]);
}

void dd_mul_dd_d ( const double *a, double b, double *c )
{
   double p1, p2;

   p1 = dd_two_prod(a[0],b,&p2);
   p2 += (a[1] * b);
   c[0] = dd_quick_two_sum(p1,p2,&c[1]);
}

void dd_mul_d_d ( double a, double b, double *c )
{
   c[0] = dd_two_prod(a,b,&c[1]);
}

void dd_mlt ( double *a, const double *b )
{
   double p1, p2;

   p1 = dd_two_prod(a[0],b[0],&p2);
   p2 += b[1] * a[0];
   p2 += b[0] * a[1];
   a[0] = dd_quick_two_sum(p1,p2,&a[1]);
}

void dd_mlt_d ( double *a, double b )
{
   double p1, p2;

   p1 = dd_two_prod(a[0],b,&p2);
   p2 += a[1] * b;
   a[0] = dd_quick_two_sum(p1,p2,&a[1]);
}

void dd_mul_pwr2 ( const double *a, double b, double *c )
{
   c[0] = a[0]*b;
   c[1] = a[1]*b;
}

/*************************** divisions ***************************/

void dd_div ( const double *a, const double *b, double *c )
{
   /* this is the accurate_div in dd_inline.h */

   double q1, q2, q3;
   double acc[2];

   q1 = a[0] / b[0];        /* approximate quotient */
   dd_mul_d_dd(q1,b,acc);   /* acc = q1*b */
   dd_sub(a,acc,c);         /* c = a - q1 * b; */
   q2 = c[0] / b[0];
   dd_mul_d_dd(q2,b,acc);   /* acc = q2*b */ 
   dd_dec(c,acc);           /* c -= (q2 * b); */
   q3 = c[0] / b[0];
   c[0] = dd_quick_two_sum(q1,q2,&c[1]);
   dd_inc_d(c,q3);          /* c = dd_real(q1, q2) + q3; */
}

void dd_div_d_dd ( double a, const double *b, double *c )
{
   double aa[2];

   aa[0] = a;
   aa[1] = 0.0;
   dd_div(aa,b,c);
}

void dd_div_dd_d ( const double *a, double b, double *c )
{
   double q1, q2, p1, p2, s, e;

   q1 = a[0] / b;                        /* approximate quotient. */
   p1 = dd_two_prod(q1,b,&p2);           /* compute a - q1 * d */
   s = dd_two_diff(a[0],p1,&e);
   e += a[1];
   e -= p2;
   q2 = (s + e) / b;                     /* get next approximation. */
   c[0] = dd_quick_two_sum(q1,q2,&c[1]); /* renormalize */
}

void dd_div_d_d ( double a, double b, double *c )
{
   double q1, q2, p1, p2, s, e;

   q1 = a / b;
   p1 = dd_two_prod(q1,b,&p2);    /* compute a - q1 * b */
   s = dd_two_diff(a,p1,&e);
   e -= p2;
   q2 = (s + e) / b;              /* get next approximation */

   c[0] = dd_quick_two_sum(q1,q2,&c[1]);
}

/************************* comparisons *******************************/

int dd_is_zero ( const double *a )
{
   return (a[0] == 0.0) ? 1 : 0;
}

int dd_is_one ( const double *a )
{
   return ((a[0] == 1.0) && (a[1] == 0.0)) ? 1 : 0;
}

int dd_is_positive ( const double *a )
{
   return (a[0] > 0.0) ? 1 : 0;
}

int dd_is_negative ( const double *a )
{
   return (a[0] < 0.0) ? 1 : 0;
}

int dd_eq_dd_d ( const double *a, double b )
{
   return ((a[0] == b) && (a[1] == 0.0)) ? 1 : 0;
}

int dd_eq ( const double *a, const double *b )
{
   return ((a[0] == b[0]) && (a[1] == b[1])) ? 1 : 0;
}

int dd_eq_d_dd ( double a, const double *b )
{
   return ((a == b[0]) && (b[1] == 0.0)) ? 1 : 0;
}

int dd_neq_dd_d ( const double *a, double b )
{
   return ((a[0] != b) || (a[1] != 0.0)) ? 1 : 0;
}

int dd_neq ( const double *a, const double *b )
{
   return ((a[0] != b[0]) || (a[1] != b[1])) ? 1 : 0;
}

int dd_neq_d_dd ( double a, const double *b )
{
   return ((a != b[0]) || (b[1] != 0.0)) ? 1 : 0;
}

int dd_gt_dd_d ( const double *a, double b )
{
   return (a[0] > b || (a[0] == b && a[1] > 0.0)) ? 1 : 0;
}

int dd_gt ( const double *a, double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] > b[1])) ? 1 : 0;
}

int dd_gt_d_dd ( double a, const double *b )
{
   return (a > b[0] || (a == b[0] && b[1] < 0.0)) ? 1 : 0;
}

int dd_geq_dd_d ( const double *a, double b )
{
   return (a[0] > b || (a[0] == b && a[1] >= 0.0)) ? 1 : 0;
}

int dd_geq ( const double *a, const double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] >= b[1])) ? 1 : 0;
}

int dd_geq_d_dd ( double a, const double *b )
{
   return (a > b[0] || (a == b[0] && b[1] <= 0.0)) ? 1 : 0;
}

int dd_lt_dd_d ( const double *a, double b )
{
   return (a[0] < b || (a[0] == b && a[1] < 0.0)) ? 1 : 0;
}

int dd_lt ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])) ? 1 : 0;
}

int dd_lt_d_dd ( double a, const double *b )
{
   return (a < b[0] || (a == b[0] && b[1] > 0.0)) ? 1 : 0;
}

int dd_leq_dd_d ( const double *a, double b )
{
   return (a[0] < b || (a[0] == b && a[1] <= 0.0)) ? 1 : 0;
}

int dd_leq ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] <= b[1])) ? 1 : 0;
}

int dd_leq_d_dd ( double a, const double *b )
{
   return (a < b[0] || (a == b[0] && b[1] >= 0.0)) ? 1 : 0;
}

/******************** squaring and power ***************************/

void dd_sqr ( const double *a, double *b )
{
   double p1, p2;

   p1 = dd_two_sqr(a[0],&p2);
   p2 += 2.0 * a[0] * a[1];
   p2 += a[1] * a[1];
   b[0] = dd_quick_two_sum(p1,p2,&b[1]);
}

void dd_sqr_d ( double a, double *b )
{
   b[0] = dd_two_sqr(a,&b[1]);
}

void dd_npwr ( const double *a, int n, double *b )
{
   if(n == 0)   /* the original crashes at 0^0 = 1 */
   {
      b[0] = 1.0; b[1] = 0.0;
   }
   else
   {
      int N = (n > 0) ? n : -n;  /* N = abs(n) */
      double s[2];

      b[0] = a[0]; s[0] = 1.0;
      b[1] = a[1]; s[1] = 0.0;

      if(N > 1)     /* use binary exponentiation */
      {
         while(N > 0)
         {
            if(N % 2 == 1) dd_mlt(s,b);
            N /= 2;
            if(N > 0) dd_sqr(b,b);
         }
      }
      else
      {
         s[0] = b[0]; s[1] = b[1];
      }
      if(n < 0) 
         dd_div_d_dd(1.0,s,b);      /* compute reciprocal */
      else
      {
         b[0] = s[0]; b[1] = s[1];
      }
   }
}

void dd_ldexp ( double *a, int n )
{
   a[0] = ldexp(a[0],n);
   a[1] = ldexp(a[1],n);
}

/*********************** exp and log ******************************/

void dd_exp ( const double *a, double *b )
{
  /* Strategy:  We first reduce the size of x by noting that

          exp(kr + m * log(2)) = 2^m * exp(r)^k

     where m and k are integers.  By choosing m appropriately
     we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
     evaluated using the familiar Taylor series.  Reducing the
     argument substantially speeds up the convergence.       */

   const double k = 512.0;
   const double inv_k = 1.0 / k;
   const double exp1hi = 2.718281828459045091e+00; /* exp(1) high word */
   const double exp1lo = 1.445646891729250158e-16; /* exp(1) low word */
   const double log2hi = 6.931471805599452862e-01; /* log(2) high word */
   const double log2lo = 2.319046813846299558e-17; /* log(2) low word */
   const double dd_eps = 4.93038065763132e-32;     /* 2^-104 */
   double m,tol;
   double r[2];
   double s[2];
   double t[2];
   double p[2];
   int i;
   double *inv_fac;
   double i_fac[30];  /* inverse factorials for Taylor expansion */
   i_fac[0] = 1.66666666666666657e-01;  i_fac[1] =  9.25185853854297066e-18;
   i_fac[2] = 4.16666666666666644e-02;  i_fac[3] =  2.31296463463574266e-18;
   i_fac[4] = 8.33333333333333322e-03;  i_fac[5] =  1.15648231731787138e-19;
   i_fac[6] = 1.38888888888888894e-03;  i_fac[7] = -5.30054395437357706e-20;
   i_fac[8] = 1.98412698412698413e-04;  i_fac[9] =  1.72095582934207053e-22;
   i_fac[10] = 2.48015873015873016e-05; i_fac[11] =  2.15119478667758816e-23;
   i_fac[12] = 2.75573192239858925e-06; i_fac[13] = -1.85839327404647208e-22;
   i_fac[14] = 2.75573192239858883e-07; i_fac[15] =  2.37677146222502973e-23;
   i_fac[16] = 2.50521083854417202e-08; i_fac[17] = -1.44881407093591197e-24;
   i_fac[18] = 2.08767569878681002e-09; i_fac[19] = -1.20734505911325997e-25;
   i_fac[20] = 1.60590438368216133e-10; i_fac[21] =  1.25852945887520981e-26;
   i_fac[22] = 1.14707455977297245e-11; i_fac[23] =  2.06555127528307454e-28;
   i_fac[24] = 7.64716373181981641e-13; i_fac[25] =  7.03872877733453001e-30;
   i_fac[26] = 4.77947733238738525e-14; i_fac[27] =  4.39920548583408126e-31;
   i_fac[28] = 2.81145725434552060e-15; i_fac[29] =  1.65088427308614326e-31;

   if(a[0] <= -709.0)
   {
      b[0] = 0.0; b[1] = 0.0; return;
   }
   if(a[0] >=  709.0)
   {
      b[0] = -1.0; b[1] = 0.0; return;
   }
   if(dd_is_zero(a) == 1)
   {
      b[0] = 1.0; b[1] = 0.0; return;
   }
   if(dd_is_one(a) == 1)
   {
      b[0] = exp1hi; b[1] = exp1lo; return;
   }
   m = floor(a[0] / log2hi + 0.5);
   /* r = mul_pwr2(a - dd_real::_log2 * m, inv_k); */
   r[0] = log2hi; r[1] = log2lo; 
   dd_mlt_d(r,m);
   dd_sub(a,r,r);
   dd_mul_pwr2(r,inv_k,r);
   dd_sqr(r,p);
   dd_mul_pwr2(p,0.5,s);   /* s = r + mul_pwr2(p, 0.5); */
   dd_inc(s,r);
   dd_mlt(p,r);            /* p *= r; */
   inv_fac = &i_fac[0];
   dd_mul(p,inv_fac,t);    /* t = p * inv_fact[0]; */
   i = 0;
   tol = inv_k*dd_eps;
   do
   {
      dd_inc(s,t);
      dd_mlt(p,r);
      inv_fac = &i_fac[2*(++i)];
      dd_mul(p,inv_fac,t);
   } while(fabs(t[0]) > tol && i < 5);

   dd_inc(s,t);            /* s += t; */
   for(i=0; i<9; i++)      /* 9 times s = mul_pwr2(s,2.0) + sqr(s); */
   {
      dd_mul_pwr2(s,2.0,p);
      dd_sqr(s,t);
      dd_add(p,t,s);
   }
   dd_add_dd_d(s,1.0,b);   /* s += 1.0; but then b stores the result */
   i = (int) m;
   dd_ldexp(b,i);
}

void dd_log ( const double *a, double *b )
{
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x)
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.

     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

   double acc[2];

   if(dd_is_one(a) == 1)
   {
      b[0] = 0.0; b[1] = 0.0; return;
   }
   if(a[0] <= 0.0)
   {
      printf("(dd_log): argument is not positive");
      b[0] = -1.0; b[0] = 0.0; return;
   }
   b[0] = log(a[0]); b[1] = 0.0;   /* initial approximation */
   acc[0] = -b[0]; acc[1] = -b[1]; /* x = x + a * exp(-x) - 1.0 */
   dd_exp(acc,acc);
   dd_mlt(acc,a);
   dd_inc(b,acc);
   dd_dec_d(b,1.0);
}

void dd_log10 ( const double *a, double *b )
{
   const double log10hi =  2.302585092994045901e+00; /* log(10) hi word */
   const double log10lo = -2.170756223382249351e-16; /* log(10) lo word */
   double logten[2];

   logten[0] = log10hi;
   logten[1] = log10lo;
   dd_log(a,b);
   dd_div(b,logten,b);   /* b = log(a)/log(10) */
}
/****************** copy, type casts, abs, and floor *****************/

void dd_copy ( const double *a, double *b )
{
   b[0] = a[0];
   b[1] = a[1];
}

int dd_to_int ( const double *a )
{
   return ((int) a[0]);
}

double dd_to_double ( const double *a )
{
   return a[0];
}

void dd_abs ( const double *a, double *b )
{
   if(a[0] < 0.0)
   {
      b[0] = -a[0];
      b[1] = -a[1];
   }
   else
   {
      b[0] = a[0];
      b[1] = a[1];
   }
}

void dd_floor ( const double *a, double *b )
{
   double hi = floor(a[0]);
   double lo = 0.0;

   if(hi == a[0])    /* high word is already integer, round low word */
   {
      lo = floor(a[1]);
      b[0] = dd_quick_two_sum(hi,lo,&b[1]);
   }
   else
   {
      b[0] = hi;
      b[1] = lo;
   }
}

void dd_nint ( const double *a, double *b )
{
   double f[2];

   dd_floor(a,f);

   if(dd_eq(a,f) == 1) 
   {
      dd_copy(a,b);
   }
   else
   {
      double c[2];

      dd_add_dd_d(a,0.5,c);
      dd_floor(c,b);
   }
}

/************************ sqrt, sin and cos ******************************/  

void dd_sqrt ( const double* a, double *b )
{
  /* Use Karp's trick: if x is an approximation to sqrt(a), then
       sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision. */
  
   if(dd_is_zero(a) == 1)
   {
      b[0] = 0.0; b[1] = 0.0;
   }
   else if(a[0] < 0.0)
   {
      b[0] = -1.0; b[0] = 0.0;
   }
   else
   {
      const double x = 1.0/sqrt(a[0]);
      double ax = a[0]*x;

      double ddax[2];
      double sqax[2];
      double y[2];
      
      dd_sqr_d(ax,sqax);
      dd_sub(a,sqax,y);     

      b[0] = dd_two_sum(ax,y[0]*(x*0.5),&b[1]);
   }
}

void dd_sin_taylor ( const double *a, double *b )
{
   if(dd_is_zero(a) == 1)
   {
      b[0] = 0.0; b[1] = 0.0;
   }
   else
   {
      const int n_inv_fact = 15;
      double i_fac[30];  /* inverse factorials for Taylor expansion */
      i_fac[0] = 1.66666666666666657e-01;  i_fac[1] =  9.25185853854297066e-18;
      i_fac[2] = 4.16666666666666644e-02;  i_fac[3] =  2.31296463463574266e-18;
      i_fac[4] = 8.33333333333333322e-03;  i_fac[5] =  1.15648231731787138e-19;
      i_fac[6] = 1.38888888888888894e-03;  i_fac[7] = -5.30054395437357706e-20;
      i_fac[8] = 1.98412698412698413e-04;  i_fac[9] =  1.72095582934207053e-22;
      i_fac[10] = 2.48015873015873016e-05; i_fac[11] =  2.15119478667758816e-23;
      i_fac[12] = 2.75573192239858925e-06; i_fac[13] = -1.85839327404647208e-22;
      i_fac[14] = 2.75573192239858883e-07; i_fac[15] =  2.37677146222502973e-23;
      i_fac[16] = 2.50521083854417202e-08; i_fac[17] = -1.44881407093591197e-24;
      i_fac[18] = 2.08767569878681002e-09; i_fac[19] = -1.20734505911325997e-25;
      i_fac[20] = 1.60590438368216133e-10; i_fac[21] =  1.25852945887520981e-26;
      i_fac[22] = 1.14707455977297245e-11; i_fac[23] =  2.06555127528307454e-28;
      i_fac[24] = 7.64716373181981641e-13; i_fac[25] =  7.03872877733453001e-30;
      i_fac[26] = 4.77947733238738525e-14; i_fac[27] =  4.39920548583408126e-31;
      i_fac[28] = 2.81145725434552060e-15; i_fac[29] =  1.65088427308614326e-31;
      const double dd_eps = 4.93038065763132e-32;     /* 2^-104 */
      const double thresh = 0.5*a[0]*dd_eps;
      double x[2], r[2], s[2], t[2], inv_fac[2];
      int i = 0;
      dd_sqr(a,x);     /* x = -sqr(a) */
      dd_minus(x);
      dd_copy(a,s);
      dd_copy(a,r);
      do
      {
         dd_mlt(r,x);             /* r *= x */
         inv_fac[0] = i_fac[i];   /* even ones are high parts */
         inv_fac[1] = i_fac[i+1]; /* odd ones are low parts */
         dd_mul(r,inv_fac,t);     /* t = r * inv_fact[i] */
         dd_inc(s,t);             /* s += t */
         i += 4;                  /* only take even terms */
      }
      while((i < 2*n_inv_fact) && (fabs(t[0]) > thresh));
      dd_copy(s,b);
   }
}

void dd_cos_taylor ( const double *a, double *b )
{
   if(dd_is_zero(a) == 1)
   {
      b[0] = 1.0; b[1] = 0.0;
   }
   else
   {
      const int n_inv_fact = 15;
      double i_fac[30];  /* inverse factorials for Taylor expansion */
      i_fac[0] = 1.66666666666666657e-01;  i_fac[1] =  9.25185853854297066e-18;
      i_fac[2] = 4.16666666666666644e-02;  i_fac[3] =  2.31296463463574266e-18;
      i_fac[4] = 8.33333333333333322e-03;  i_fac[5] =  1.15648231731787138e-19;
      i_fac[6] = 1.38888888888888894e-03;  i_fac[7] = -5.30054395437357706e-20;
      i_fac[8] = 1.98412698412698413e-04;  i_fac[9] =  1.72095582934207053e-22;
      i_fac[10] = 2.48015873015873016e-05; i_fac[11] =  2.15119478667758816e-23;
      i_fac[12] = 2.75573192239858925e-06; i_fac[13] = -1.85839327404647208e-22;
      i_fac[14] = 2.75573192239858883e-07; i_fac[15] =  2.37677146222502973e-23;
      i_fac[16] = 2.50521083854417202e-08; i_fac[17] = -1.44881407093591197e-24;
      i_fac[18] = 2.08767569878681002e-09; i_fac[19] = -1.20734505911325997e-25;
      i_fac[20] = 1.60590438368216133e-10; i_fac[21] =  1.25852945887520981e-26;
      i_fac[22] = 1.14707455977297245e-11; i_fac[23] =  2.06555127528307454e-28;
      i_fac[24] = 7.64716373181981641e-13; i_fac[25] =  7.03872877733453001e-30;
      i_fac[26] = 4.77947733238738525e-14; i_fac[27] =  4.39920548583408126e-31;
      i_fac[28] = 2.81145725434552060e-15; i_fac[29] =  1.65088427308614326e-31;
      const double dd_eps = 4.93038065763132e-32;     /* 2^-104 */
      const double thresh = 0.5*dd_eps;
      double x[2], r[2], s[2], t[2], inv_fac[2];
      int i = 1;
      dd_sqr(a,x);                /* x = -sqr(a) */
      dd_minus(x);
      dd_copy(x,r);
      dd_mul_pwr2(r, 0.5, s);     /* s = 1.0 + mul_pwr2(r, 0.5) */
      dd_inc_d(s, 1.0);
      do
      {
         dd_mlt(r,x);             /* r *= x */
         inv_fac[0] = i_fac[i+1]; /* even ones are high parts */
         inv_fac[1] = i_fac[i+2]; /* odd ones are low parts */
         dd_mul(r, inv_fac, t);   /* t = r * inv_fact[i] */
         dd_inc(s,t);             /* s += t */
         i += 4;                  /* only take the odd terms */
      }
      while((i < 2*n_inv_fact) && (fabs(t[0]) > thresh));
      dd_copy(s,b);
   }
}

void dd_sincos_taylor ( const double *a, double *sin_a, double *cos_a )
{
   if(dd_is_zero(a) == 1)
   {
      sin_a[0] = 0.0; sin_a[1] = 0.0;
      cos_a[0] = 1.0; cos_a[1] = 0.0;
   }
   else
   {
      double tmp[2];

      dd_sin_taylor(a, sin_a);
      dd_sqr(sin_a, tmp);       /* tmp = sqr(sin_a) */
      dd_minus(tmp);            /* tmp = -sqr(sin_a) */
      dd_inc_d(tmp, 1.0);       /* tmp = 1.0 - sqr(sin_a) */
      dd_sqrt(tmp, cos_a);      /* cos_a = sqrt(1.0 - sqr(sin_a) */
   }
}

void dd_reduce_modulo_2pi
 ( const double *x, double *t, int *j, int *k, int *abs_k, int* fail )
{
   double q,y[2],z[2],r[2],twopi[2],pi2[2],pi16[2];

   const double twopi_hi = 6.283185307179586232e+00; // high part of 2*pi
   const double twopi_lo = 2.449293598294706414e-16; // low part of 2*pi
   const double pi2_hi = 1.570796326794896558e+00; // high part of pi/2
   const double pi2_lo = 6.123233995736766036e-17; // low part of pi/2
   const double pi16_hi = 1.963495408493620697e-01; // high part of pi/16
   const double pi16_lo = 7.654042494670957545e-18; // low part of pi/16

   *j = 0; *k = 0; *abs_k = 0; *fail = 0;

   twopi[0] = twopi_hi; twopi[1] = twopi_lo;
   pi2[0] = pi2_hi;     pi2[1] = pi2_lo;
   pi16[0] = pi16_hi;   pi16[1] = pi16_lo;

   dd_div(x,twopi,y);

   dd_nint(y,z);
   dd_mul(twopi,z,y);
   dd_sub(x,y,r);

   q = floor(r[0]/pi2_hi+0.5);
   dd_mul_dd_d(pi2,q,y);
   dd_sub(r,y,t);
   *j = (int) q;

   if((*j < -2 ) || (*j > 2))
   {
      printf("dd_sin: cannot reduce module pi/2");
      *fail = 1;
   }
   else
   {
      q = floor(t[0]/pi16_hi + 0.5);
      dd_mul_dd_d(pi16,q,y);
      dd_dec(t,y);
      *k = (int) q;
      if(*k < 0)
         *abs_k = -(*k);
      else
         *abs_k = *k;
      if(*abs_k > 4)
      {
         printf("dd_sin: cannot reduce module pi/16");
         *fail = 1;
      }
      else
         *fail = 0;
   }
}

void sincostables
 ( double *sintabhi, double *sintablo, double *costabhi, double *costablo )
/*
 * Defines tables of sin(k*pi/16) and cos(k*pi/16),
 * in high parts ending with hi and low parts ending with lo.
 * All arrays have space for four doubles. */
{
   const double sin_t0_hi = 1.950903220161282758e-01;
   const double sin_t0_lo = -7.991079068461731263e-18;
   const double sin_t1_hi = 3.826834323650897818e-01;
   const double sin_t1_lo = -1.005077269646158761e-17;
   const double sin_t2_hi = 5.555702330196021776e-01;
   const double sin_t2_lo = 4.709410940561676821e-17;
   const double sin_t3_hi = 7.071067811865475727e-01;
   const double sin_t3_lo = -4.833646656726456726e-17;
   const double cos_t0_hi = 9.807852804032304306e-01;
   const double cos_t0_lo = 1.854693999782500573e-17;
   const double cos_t1_hi = 9.238795325112867385e-01;
   const double cos_t1_lo = 1.764504708433667706e-17;
   const double cos_t2_hi = 8.314696123025452357e-01;
   const double cos_t2_lo = 1.407385698472802389e-18;
   const double cos_t3_hi = 7.071067811865475727e-01;
   const double cos_t3_lo = -4.833646656726456726e-17;

   sintabhi[0] = sin_t0_hi; sintablo[0] = sin_t0_lo;
   sintabhi[1] = sin_t1_hi; sintablo[1] = sin_t1_lo;
   sintabhi[2] = sin_t2_hi; sintablo[2] = sin_t2_lo;
   sintabhi[3] = sin_t3_hi; sintablo[3] = sin_t3_lo;
   costabhi[0] = cos_t0_hi; costablo[0] = cos_t0_lo;
   costabhi[1] = cos_t1_hi; costablo[1] = cos_t1_lo;
   costabhi[2] = cos_t2_hi; costablo[2] = cos_t2_lo;
   costabhi[3] = cos_t3_hi; costablo[3] = cos_t3_lo;
}

void dd_sin ( const double *a, double *sin_a )
{
   if(dd_is_zero(a) == 1)
   {
      sin_a[0] = 0.0; sin_a[1] = 0.0;
   }
   else
   {
      double t[2];
      int j,k,abs_k,fail;

      dd_reduce_modulo_2pi(a,t,&j,&k,&abs_k,&fail);
      if(fail == 1)
      {
         sin_a[0] = -2.0; sin_a[1] = 0.0;
      }
      else if(k == 0)
      {
         if(j == 0)
         {
            dd_sin_taylor(t,sin_a);
         }
         else if(j == 1)
         {
            dd_cos_taylor(t,sin_a);
         }
         else if(j == -1)
         {
            dd_cos_taylor(t,sin_a);
            dd_minus(sin_a);
         }
         else
         {
            dd_sin_taylor(t,sin_a);
            dd_minus(sin_a);
         }
      }
      else
      {
         double u[2],v[2],w[2],t_sin[2],t_cos[2];
         double sintabhi[4],sintablo[4];
         double costabhi[4],costablo[4];

         sincostables(sintabhi,sintablo,costabhi,costablo);

         u[0] = costabhi[abs_k-1]; u[1] = costablo[abs_k-1];
         v[0] = sintabhi[abs_k-1]; v[1] = sintablo[abs_k-1];
         dd_sincos_taylor(t,t_sin,t_cos);
         if(j == 0)
         {
            dd_mul(u,t_sin,sin_a);
            dd_mul(v,t_cos,w);
            if(k > 0)
               dd_inc(sin_a,w);
            else
               dd_dec(sin_a,w);
         }
         else if(j == 1)
         {
            dd_mul(u,t_cos,sin_a);
            dd_mul(v,t_sin,w);
            if(k > 0)
               dd_dec(sin_a,w);
            else
               dd_inc(sin_a,w);
         }
         else if(j == -1)
         {
            dd_mul(v,t_sin,sin_a);
            dd_mul(u,t_cos,w);
            if(k > 0)
            {
               dd_dec(sin_a,w);
            }
            else
            { 
               dd_minus(sin_a);
               dd_dec(sin_a,w);
            }
         }
         else
         {
            dd_mul(v,t_cos,sin_a);
            dd_mul(u,t_sin,w);
            if(k > 0)
            {
               dd_minus(sin_a);
               dd_dec(sin_a,w);
            }
            else
            {
               dd_dec(sin_a,w);
            }
         }
      }
   }
}

void dd_cos ( const double *a, double *cos_a )
{
   if(dd_is_zero(a) == 1)
   {
      cos_a[0] = 1.0; cos_a[1] = 0.0;
   }
   else
   {
      double t[2];
      int j,k,abs_k,fail;
 
      dd_reduce_modulo_2pi(a,t,&j,&k,&abs_k,&fail);
      if(fail == 1)
      {
         cos_a[0] = -2.0; cos_a[1] = 0.0;
      }
      else if(k == 0)
      {
         if(j == 0)
         {
            dd_cos_taylor(t,cos_a);
         }
         else if(j == 1)
         {
            dd_sin_taylor(t,cos_a);
            dd_minus(cos_a);
         }
         else if(j == -1)
         {
            dd_sin_taylor(t,cos_a);
         }
         else
         {
            dd_cos_taylor(t,cos_a);
            dd_minus(cos_a);
         }
      }
      else
      {
         double u[2],v[2],w[2],t_sin[2],t_cos[2];
         double sintabhi[4],sintablo[4];
         double costabhi[4],costablo[4];

         sincostables(sintabhi,sintablo,costabhi,costablo);

         u[0] = costabhi[abs_k-1]; u[1] = costablo[abs_k-1];
         v[0] = sintabhi[abs_k-1]; v[1] = sintablo[abs_k-1];
         dd_sincos_taylor(t,t_sin,t_cos);
         if(j == 0)
         {
            dd_mul(u,t_cos,cos_a);
            dd_mul(v,t_sin,w);
            if(k > 0)
               dd_dec(cos_a,w);
            else
               dd_inc(cos_a,w);
         }
         else if(j == 1)
         {
            dd_mul(v,t_cos,cos_a);
            dd_mul(u,t_sin,w);
            if(k > 0)
            {
               dd_minus(cos_a);
               dd_dec(cos_a,w);
            }
            else
            {
               dd_dec(cos_a,w);
            }
         }
         else if(j == -1)
         {
            dd_mul(u,t_sin,cos_a);
            dd_mul(v,t_cos,w);
            if(k > 0)
               dd_inc(cos_a,w);
            else
               dd_dec(cos_a,w);
         }
         else
         {
            dd_mul(v,t_sin,cos_a);
            dd_mul(u,t_cos,w);
            if(k > 0)
            {
               dd_dec(cos_a,w);
            }
            else
            {
               dd_minus(cos_a);
               dd_dec(cos_a,w);
            }
         }
      }
   }
}

/********************* input/output operations *************************/

int dd_read ( const char *s, double *a )
{
   const char *p = s;
   char ch;
   int sign = 0;
   int point = -1;
   int nd = 0;
   int e = 0;
   int done = 0;
   int nread;

   while (*p == ' ') p++;    /* skip leading spaces */

   a[0] = 0.0; a[1] = 0.0;
   while ((done == 0) && (ch = *p) != '\0')
   {
      if (ch >= '0' && ch <= '9')
      {
         int d = ch - '0';
         double cd = (double) d;

         dd_mlt_d(a,10.0);  /* a *= 10.0 */
         dd_inc_d(a,cd);    /* a += (double) d */
         nd++;
      } 
      else
      {
         switch (ch)
         {
            case '.':
               if (point >= 0) return -1;
               point = nd; break;
            case '-':
            case '+':
               if (sign != 0 || nd > 0) return -1;
               sign = (ch == '-') ? -1 : 1; break;
            case 'E':
            case 'e':
               nread = sscanf(p+1, "%d", &e);
               done = 1;
               if (nread != 1) return -1; break;
            default:
               return -1;
         }
      }
      p++;
   }
   if (point >= 0) e -= (nd - point);

   if(e != 0)  /* r *= (dd_real(10.0) ^ e); */
   {
      double acc[2];

      acc[0] = 10.0;
      acc[1] = 0.0;
      dd_npwr(acc,e,acc);
      dd_mlt(a,acc);
   }
   if(sign == -1) dd_minus(a);  /* a = (sign == -1) ? -r : r; */

   return 0;
}

void dd_to_digits ( const double *a, char *s, int *expn, int precision )
{
   int D = precision + 1;  /* number of digits to compute */
   double acc[2];          /* for the r in the original code */
   double tmp[2];          /* to hold temporary results */
   int i,d,e;              /* e is exponent */

   dd_abs(a,acc);
   
   if(a[0] == 0.0)         /* a equals zero */
   {
      expn = 0;
      for(i = 0; i < precision; i++) s[i] = '0'; return;
   }
   e = (int) (floor(log10(fabs(a[0]))));  /* approximate exponent */
   if (e < -300)
   {
      tmp[0] = 10.0; tmp[1] = 0.0;
      dd_npwr(tmp,300,tmp);
      dd_mlt(acc,tmp);                 /* r *= dd_real(10.0) ^ 300; */
      tmp[0] = 10.0; tmp[1] = 0.0;
      dd_npwr(tmp,e+300,tmp);
      dd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ (e + 300); */
   }
   else if(e > 300)
   {
      dd_ldexp(acc,-53);               /* r = ldexp(r, -53); */
      tmp[0] = 10.0; tmp[1] = 0.0;
      dd_npwr(tmp,e,tmp);
      dd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ e; */
      dd_ldexp(acc,53);                /* r = ldexp(r, 53); */
   }
   else
   {
      tmp[0] = 10.0; tmp[1] = 0.0;
      dd_npwr(tmp,e,tmp);
      dd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ e; */
   }
   if(dd_gt_dd_d(acc,10.0) == 1)  /* fix exponent if we are off by one */
   {
      dd_div_dd_d(acc,10.0,acc);       /* r /= 10.0; */
      e++;
   }
   else if(dd_lt_dd_d(acc,1.0) == 1)   /* r < 1.0) */
   {
      dd_mlt_d(acc,10.0);              /* r *= 10.0; */
      e--;
   }
   if((dd_geq_dd_d(acc,10.0) == 1) || (dd_lt_dd_d(acc,1.0) == 1))
   {                 /* (r >= 10.0 || r < 1.0) */
      printf("dd_to_digits: cannot compute exponent"); return;
   }
   for(i=0; i<D; i++)                  /* extract the digits */
   {
      d = ((int) acc[0]);
      dd_dec_d(acc,(double) d);        /* r -= d; */
      dd_mlt_d(acc,10.0);              /* r *= 10.0; */
      s[i] = ((char) (d + '0'));
   }
   for(i=D-1; i>0; i--)    /* fix out of range digits */
   {
      if(s[i] < '0')
      {
         s[i-1]--;
         s[i] += 10;
      }
      else if(s[i] > '9')
      {
         s[i-1]++;
         s[i] -= 10;
      }
   }
   if(s[0] <= '0') 
   {
      printf("dd_to_digits: nonpositive leading digit"); return;
   }
   if(s[D-1] >= '5')   /* round, handle carry */
   {
      s[D-2]++;
      i = D-2;
      while (i > 0 && s[i] > '9')
      {
         s[i] -= 10;
         s[--i]++;
      }
   }
   if(s[0] > '9')    /* if first digit is 10, shift everything */
   {
      e++;
      for(i=precision; i>=2; i--) s[i] = s[i-1];
      s[0] = '1';
      s[1] = '0';
   }
   s[precision] = 0;
   *expn = e;
}

void dd_to_string ( const double *a, char *s, int precision, int width,
                    int fixed, int showpos, int uppercase, char fill,
                    int *endstring )
{
/* note: nan and inf are ignored ...*/

   int sgn = 1;
   int i, e = 0;
   char *p = s;
   int cnt = 0; /* counts characters written to string */
   int off,d;
   double acc[2];

   if(dd_is_negative(a) == 1)
   {
      *(p++) = '-'; cnt++;
   }
   else if(showpos == 1)
   {
      *(p++) = '+'; cnt++;
   }
   else
      sgn = 0;

   if(dd_is_zero(a) == 1)
   {
      *(p++) = '0'; cnt++;
      if(precision > 0)
      {
         *(p++) = '.'; cnt++;
         for(i=0; i<precision; i++,cnt++) *(p++) = '0';
      }
   }
   else  /* nonzero case */
   {
      if(fixed != 1)
         off = 1;
      else
      {
         dd_abs(a,acc);
         dd_log10(acc,acc);
         dd_floor(acc,acc);
         off = dd_to_int(acc) + 1;
      }
      d = precision + off;
      if((fixed == 1) && (d <= 0))
      {
         *(p++) = '0'; cnt++;
         if(precision > 0)
         {
            *(p++) = '.'; cnt++;
            for(i=0; i<precision; i++,cnt++) *(p++) = '0';
         }
      }
      else
      {
         char* t = (char*)calloc(d+1,sizeof(char)); // char t[d+1]
         int j;

         dd_to_digits(a,t,&e,d);
         if(fixed == 1)
         {
            if (off > 0)
            {
               for (i=0; i<off; i++,cnt++) *(p++) = t[i];
               if(precision > 0) 
               {
                  *(p++) = '.'; cnt++;
                  for(j=0; j<precision; j++,i++,cnt++) *(p++) = t[i];
               }
            }
            else
            {
               *(p++) = '0'; cnt++;
               *(p++) = '.'; cnt++;
               if(off < 0) 
                  for(i=0; i<(-off); i++,cnt++) *(p++) = '0';
               for(i=0; i<d; i++,cnt++) *(p++) = t[i];
            }
         } 
         else
         {
            *(p++) = t[0]; cnt++;
            if(precision > 0)
            {
               *(p++) = '.'; cnt++;
            }
            for(i=1; i<=precision; i++,cnt++) *(p++) = t[i];
         }
         free(t);
      }
   }
   if(fixed == 0)   /* fill in the exponent part */
   {
      *(p++) = (uppercase == 1) ? 'E' : 'e'; cnt++;
      *(p++) = (e < 0 ? '-' : '+'); cnt++;
      if(e<0) e = -e;
      if (e >= 100)
      {
         i = e/100;
         *(p++) = '0' + i; cnt++;
         e -= 100*i;
      }
      i = e/10;
      *(p++) = '0' + i; cnt++;
      e -= 10*i;
      *(p++) = '0' + e; cnt++;
   }
   if(cnt >= width)
   {
      *(p++) = '\0';
      *endstring = cnt + 1;
   }
   else   /* fill in the blanks */
   {
      d = width - cnt;
      for(i=0; i<cnt; i++) s[i+d] = s[i];
      for(i=0; i<d; i++) s[i] = fill;
      s[width] = '\0';
      *endstring = width;
   }
}

void dd_write ( const double *a, int precision )
{
   char* s = (char*)calloc(precision+10,sizeof(char)); // char s[precision+10];
   int s_end;

   dd_to_string(a,s,precision,0,0,0,1,' ',&s_end);
   printf("%s",s);
   free(s);
}
