/* file: quad_double.c */

/* This file contains the corresponding C code for the functions
   with prototypes declared in the quad_double.h file. */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"double_double.h"

/* Part I: basic functions from qd_inline.h */

void qd_quick_renorm ( double *c0, double *c1, double *c2,
                       double *c3, double *c4 )
{
   double s,t0,t1,t2,t3;

   s = dd_quick_two_sum(*c3,*c4,&t3);
   s = dd_quick_two_sum(*c2,s,&t2);
   s = dd_quick_two_sum(*c1,s,&t1);
   *c0 = dd_quick_two_sum(*c0,s,&t0);
   s = dd_quick_two_sum(t2,t3,&t2);
   s = dd_quick_two_sum(t1,s,&t1);
   *c1 = dd_quick_two_sum(t0,s,&t0);
   s = dd_quick_two_sum(t1,t2,&t1);
   *c2 = dd_quick_two_sum(t0,s,&t0);
   *c3 = t0 + t1;
}

void qd_renorm4 ( double *c0, double *c1, double *c2, double *c3 )
{
   double s0, s1;
   double s2 = 0.0;
   double s3 = 0.0;

   s0 = dd_quick_two_sum(*c2,*c3,c3);
   s0 = dd_quick_two_sum(*c1,s0,c2);
   *c0 = dd_quick_two_sum(*c0,s0,c1);

   s0 = *c0; s1 = *c1;
   if(s1 != 0.0)
   {
      s1 = dd_quick_two_sum(s1,*c2,&s2);
      if(s2 != 0.0)
         s2 = dd_quick_two_sum(s2,*c3,&s3);
      else
         s1 = dd_quick_two_sum(s1,*c3,&s2);
   }
   else
   {
      s0 = dd_quick_two_sum(s0,*c2,&s1);
      if(s1 != 0.0)
         s1 = dd_quick_two_sum(s1,*c3,&s2);
      else
         s0 = dd_quick_two_sum(s0,*c3,&s1);
   }
   *c0 = s0; *c1 = s1; *c2 = s2; *c3 = s3;
}

void qd_renorm5 ( double *c0, double *c1, double *c2,
                  double *c3, double *c4 )
{
   double s0, s1;
   double s2 = 0.0;
   double s3 = 0.0;

   s0 = dd_quick_two_sum(*c3,*c4,c4);
   s0 = dd_quick_two_sum(*c2,s0,c3);
   s0 = dd_quick_two_sum(*c1,s0,c2);
   *c0 = dd_quick_two_sum(*c0,s0,c1);
   s0 = *c0; s1 = *c1;
   s0 = dd_quick_two_sum(*c0,*c1,&s1);
   if(s1 != 0.0) 
   {
      s1 = dd_quick_two_sum(s1,*c2,&s2);
      if(s2 != 0.0) 
      {
         s2 = dd_quick_two_sum(s2,*c3,&s3);
         if(s3 != 0.0)
            s3 += *c4;
         else
            s2 += *c4;
      }
      else
      {
         s1 = dd_quick_two_sum(s1,*c3,&s2);
         if(s2 != 0.0)
            s2 = dd_quick_two_sum(s2,*c4,&s3);
         else
            s1 = dd_quick_two_sum(s1,*c4,&s2);
      }
   }
   else
   {
      s0 = dd_quick_two_sum(s0,*c2,&s1);
      if(s1 != 0.0)
      {
         s1 = dd_quick_two_sum(s1,*c3,&s2);
         if(s2 != 0.0)
            s2 = dd_quick_two_sum(s2,*c4,&s3);
         else
            s1 = dd_quick_two_sum(s1,*c4,&s2);
      }
      else
      {
         s0 = dd_quick_two_sum(s0,*c3,&s1);
         if(s1 != 0.0)
            s1 = dd_quick_two_sum(s1,*c4,&s2);
         else
            s0 = dd_quick_two_sum(s0,*c4,&s1);
      }
   }
   *c0 = s0; *c1 = s1; *c2 = s2; *c3 = s3;
}

/************************* additions *****************************/

void qd_three_sum ( double *a, double *b, double *c )
{
   double t1, t2, t3;

   t1 = dd_two_sum(*a,*b,&t2);
   *a = dd_two_sum(*c,t1,&t3);
   *b = dd_two_sum(t2,t3,c);
}

void qd_three_sum2 ( double *a, double *b, double *c )
{
   double t1, t2, t3;

   t1 = dd_two_sum(*a,*b,&t2);
   *a = dd_two_sum(*c,t1,&t3);
   *b = t2 + t3;
}

void qd_add_d ( const double *a, double b, double *c )
{
   double c0, c1, c2, c3, e;

   c[0] = dd_two_sum(a[0],b,&e);
   c[1] = dd_two_sum(a[1],e,&e);
   c[2] = dd_two_sum(a[2],e,&e);
   c[3] = dd_two_sum(a[3],e,&e);
   qd_renorm5(&c[0],&c[1],&c[2],&c[3],&e);
}

void qd_add_dd ( const double *a, const double *b, double *c )
{
   double s0, s1, s2, s3, t0, t1;

   s0 = dd_two_sum(a[0],b[0],&t0);
   s1 = dd_two_sum(a[1],b[1],&t1);
   s1 = dd_two_sum(s1,t0,&t0);
   s2 = a[2];
   qd_three_sum(&s2,&t0,&t1);
   s3 = dd_two_sum(t0,a[3],&t0);
   t0 += t1;
   qd_renorm5(&s0,&s1,&s2,&s3,&t0);
   c[0] = s0; c[1] = s1;
   c[2] = s2; c[3] = s3;
}

double qd_quick_three_accum ( double *a, double *b, double c)
{
   double s;
   int za, zb;

   s = dd_two_sum(*b,c,b);
   s = dd_two_sum(*a,s,a);

   za = (*a != 0.0) ? 1 : 0;
   zb = (*b != 0.0) ? 1 : 0;

   if((za == 1) && (zb == 1)) return s;

   if (zb == 0)
   {
      *b = *a; *a = s;
   }
   else
      *a = s;

   return 0.0;
}

void qd_add ( const double *a, const double *b, double *c )
{
/* note: this is the accurate version satifying IEEE error bound */

   int i, j, k;
   double s, t;
   double u, v;   /* double-length accumulator */
   double r[4];

   r[0] = 0.0; r[1] = 0.0; r[2] = 0.0; r[3] = 0.0;
   i = 0; j = 0; k = 0;
   if(fabs(a[i]) > fabs(b[j]))
      u = a[i++];
   else
      u = b[j++];
   if(fabs(a[i]) > fabs(b[j]))
      v = a[i++];
   else
      v = b[j++];
   u = dd_quick_two_sum(u,v,&v);
   while (k < 4)
   {
      if(i >= 4 && j >= 4)
      {
         r[k] = u;
         if (k < 3) r[++k] = v;
         break;
      }
      if(i >= 4)
         t = b[j++];
      else if(j >= 4)
         t = a[i++];
      else if (fabs(a[i]) > fabs(b[j])) 
         t = a[i++];
      else
         t = b[j++];
      s = qd_quick_three_accum(&u,&v,t);
      if(s != 0.0) r[k++] = s;
   }
   for(k=i; k<4; k++) r[3] += a[k];    /* add the rest */
   for(k=j; k<4; k++) r[3] += b[k];
   qd_renorm4(&r[0],&r[1],&r[2],&r[3]);
   c[0] = r[0]; c[1] = r[1]; c[2] = r[2]; c[3] = r[3];
}

/****** constructor, copy, abs, type casts, and unary minus ***********/

void qd_real ( double a, double *b )
{
   b[0] = a; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0;
}

void qd_copy ( const double *a, double *b )
{
   b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; b[3] = a[3];
}

void qd_abs ( const double *a, double *b )
{
   if(a[0] < 0.0)
   {
      b[0] = -a[0]; b[1] = -a[1]; b[2] = -a[2]; b[3] = -a[3];
   }
   else
   {
      b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; b[3] = a[3];
   }
}

int qd_to_int ( const double *a )
{
   return ((int) a[0]);
}

double qd_to_double ( const double *a )
{
   return a[0];
}

void qd_floor ( const double *a, double *b )
{
   double x0, x1, x2, x3;
  
   x1 = 0.0; x2 = 0.0; x3 = 0.0;
   x0 = floor(a[0]);

   if(x0 == a[0])
   {
      x1 = floor(a[1]);
      if(x1 == a[1])
      {
         x2 = floor(a[2]);
         if(x2 == a[2]) x3 = floor(a[3]);
      }
      qd_renorm4(&x0,&x1,&x2,&x3);
      b[0] = x0; b[1] = x1; b[2] = x2; b[3] = x3;
   }
   else
      b[0] = x0; b[1] = x1; b[2] = x2; b[3] = x3;
}

void qd_minus ( double *a )
{
   a[0] = -a[0]; a[1] = -a[1]; a[2] = -a[2]; a[3] = -a[3];
}

/*********************** subtractions **************************/

void qd_sub ( const double *a, const double *b, double *c )
{
   double mb[4];

   qd_copy(b,mb);
   qd_minus(mb);
   qd_add(a,mb,c);
}

void qd_sub_qd_d ( const double *a, double b, double *c )
{
   qd_add_d(a,-b,c);
}

void qd_sub_d_qd ( double a, const double *b, double *c )
{
   double mb[4];

   qd_copy(b,mb);
   qd_minus(mb);
   qd_add_d(mb,a,c);
}

void qd_sub_qd_dd ( const double *a, const double *b, double *c )
{
   double mb[2];

   dd_copy(b,mb);
   dd_minus(mb);
   qd_add_dd(a,mb,c);
}

void qd_sub_dd_qd ( const double *a, const double *b, double *c )
{
   double mb[4];

   qd_copy(b,mb);
   qd_minus(mb);
   qd_add_dd(mb,a,c);
}

/************************* multiplications *************************/

void qd_mul_pwr2 ( const double *a, double b, double *c )
{
   c[0] = a[0]*b; c[1] = a[1]*b; c[2] = a[2]*b; c[3] = a[3]*b;  
}

void qd_mul_d ( const double *a, double b, double *c )
{
   double p0, p1, p2, p3, q0, q1, q2, s4;

   p0 = dd_two_prod(a[0],b,&q0); 
   p1 = dd_two_prod(a[1],b,&q1); 
   p2 = dd_two_prod(a[2],b,&q2); 
   p3 = a[3] * b; 
   c[0] = p0;
   c[1] = dd_two_sum(q0,p1,&c[2]); 
   qd_three_sum(&c[2],&q1,&p2); 
   qd_three_sum2(&q1,&q2,&p3); 
   c[3] = q1;
   s4 = q2 + p2;
   qd_renorm5(&c[0],&c[1],&c[2],&c[3],&s4); 
}

void qd_mul_dd ( const double *a, const double *b, double *c )
{
/* a0 * b0                        0    
        a0 * b1                   1    
        a1 * b0                   2    
             a1 * b1              3    
             a2 * b0              4    
                  a2 * b1         5    
                  a3 * b0         6    
                       a3 * b1    7 
*/ 
   double p0, p1, p2, p3, p4;
   double q0, q1, q2, q3, q4;
   double s0, s1, s2;
   double t0, t1;

   p0 = dd_two_prod(a[0],b[0],&q0); 
   p1 = dd_two_prod(a[0],b[1],&q1); 
   p2 = dd_two_prod(a[1],b[0],&q2); 
   p3 = dd_two_prod(a[1],b[1],&q3); 
   p4 = dd_two_prod(a[2],b[0],&q4); 
   qd_three_sum(&p1,&p2,&q0); 
   qd_three_sum(&p2,&p3,&p4);       /* five-three-sum */
   q1 = dd_two_sum(q1,q2,&q2); 
   s0 = dd_two_sum(p2,q1,&t0); 
   s1 = dd_two_sum(p3,q2,&t1); 
   s1 = dd_two_sum(s1,t0,&t0); 
   s2 = t0 + t1 + p4;
   p2 = s0;
   p3 = a[2] * b[0] + a[3] * b[1] + q3 + q4;
   qd_three_sum2(&p3,&q0,&s1);
   p4 = q0 + s2;
   qd_renorm5(&p0,&p1,&p2,&p3,&p4);
   c[0] = p0; c[1] = p1; c[2] = p2; c[3] = p3;
}

void qd_mul ( const double *a, const double *b, double *c )
{
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5
                  a0 * b3     6
                  a1 * b2     7
                  a2 * b1     8
                  a3 * b0     9 
note: accurate version satisfying IEEE error bound
*/
   double p0, p1, p2, p3, p4, p5;
   double q0, q1, q2, q3, q4, q5;
   double p6, p7, p8, p9;
   double q6, q7, q8, q9;
   double r0, r1;
   double t0, t1;
   double s0, s1, s2;

   p0 = dd_two_prod(a[0],b[0],&q0);
   p1 = dd_two_prod(a[0],b[1],&q1);
   p2 = dd_two_prod(a[1],b[0],&q2);
   p3 = dd_two_prod(a[0],b[2],&q3);
   p4 = dd_two_prod(a[1],b[1],&q4);
   p5 = dd_two_prod(a[2],b[0],&q5);
   qd_three_sum(&p1,&p2,&q0);          /* start accumulation */
   qd_three_sum(&p2,&q1,&q2);          /* six-three sum */
   qd_three_sum(&p3,&p4,&p5);          /* of p2, q1, q2, p3, p4, p5 */
   s0 = dd_two_sum(p2,p3,&t0);         /* compute (s0, s1, s2) */
   s1 = dd_two_sum(q1,p4,&t1);         /*  = (p2, q1, q2) + (p3, p4, p5) */
   s2 = q2 + p5;
   s1 = dd_two_sum(s1,t0,&t0);
   s2 += (t0 + t1);
   p6 = dd_two_prod(a[0],b[3],&q6);    /* O(eps^3) order terms */
   p7 = dd_two_prod(a[1],b[2],&q7);
   p8 = dd_two_prod(a[2],b[1],&q8);
   p9 = dd_two_prod(a[3],b[0],&q9);
   q0 = dd_two_sum(q0,q3,&q3);         /* nine-two-sum of q0, s1, q3, */
   q4 = dd_two_sum(q4,q5,&q5);         /* q4, q5, p6, p7, p8, p9 */
   p6 = dd_two_sum(p6,p7,&p7);
   p8 = dd_two_sum(p8,p9,&p9);
   t0 = dd_two_sum(q0,q4,&t1);         /* compute (t0, t1) */
   t1 += (q3 + q5);                    /* = (q0, q3) + (q4, q5) */
   r0 = dd_two_sum(p6,p8,&r1);         /* compute (r0, r1) */
   r1 += (p7 + p9);                    /* = (p6, p7) + (p8, p9) */
   q3 = dd_two_sum(t0,r0,&q4);         /* compute (q3, q4) */
   q4 += (t1 + r1);                    /* = (t0, t1) + (r0, r1) */
   t0 = dd_two_sum(q3,s1,&t1);         /* compute (t0, t1) */
   t1 += q4;                           /* = (q3, q4) + s1 */
                                       /* O(eps^4) terms -- nine-one-sum */
   t1 += a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + q6 + q7 + q8 + q9 + s2;
   qd_renorm5(&p0,&p1,&s0,&t0,&t1);
   c[0] = p0; c[1] = p1; c[2] = s0; c[3] = t0;
}

/********************** comparisons *****************************/

int qd_is_zero ( const double *a )
{
   return (a[0] == 0.0) ? 1 : 0;
}

int qd_is_one ( const double *a )
{
   return ((a[0] == 1.0) && (a[1] == 0.0) && (a[2] == 0.0)
                         && (a[3] == 0.0)) ? 1 : 0;
}

int qd_is_positive ( const double *a )
{
   return (a[0] > 0.0) ? 1 : 0;
}

int qd_is_negative ( const double *a )
{
   return (a[0] < 0.0) ? 1 : 0;
}

int qd_eq_qd_d ( const double *a, double b )
{
    return ((a[0] == b) && (a[1] == 0.0)
         && (a[2] == 0.0) && (a[3] == 0.0)) ? 1 : 0;
}

int qd_eq_qd_dd ( const double *a, const double *b )
{
    return ((a[0] == b[0]) && (a[1] == b[1])
         && (a[2] == 0.0) && (a[3] == 0.0)) ? 1 : 0;
}

int qd_eq ( const double *a, const double *b )
{
    return ((a[0] == b[0]) && (a[1] == b[1])
         && (a[2] == b[2]) && (a[3] == b[3])) ? 1 : 0;
}

int qd_lt_qd_d ( const double *a, double b )
{
   return (a[0] < b || (a[0] == b && a[1] < 0.0)) ? 1 : 0;
}

int qd_lt_qd_dd ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
                       || (a[1] == b[1] && a[2] < 0.0)) ? 1 : 0;
}

int qd_lt ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
                       || (a[1] == b[1] && a[2] < b[2])
                       || (a[2] == b[2] && a[3] < b[3])) ? 1 : 0;
}

int qd_leq_qd_d ( const double *a, double b )
{
   return (a[0] < b || (a[0] == b && a[1] <= 0.0)) ? 1 : 0;
}

int qd_leq_qd_dd ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
                       || (a[1] == b[1] && a[2] <= 0.0)) ? 1 : 0;
}

int qd_leq ( const double *a, const double *b )
{
   return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
                       || (a[1] == b[1] && a[2] < b[2])
                       || (a[2] == b[2] && a[3] <= b[3])) ? 1 : 0;
}

int qd_gt_qd_d ( const double *a, double b )
{
   return (a[0] > b || (a[0] == b && a[1] > 0.0)) ? 1 : 0;
}

int qd_gt_qd_dd ( const double *a, const double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] > b[1])
                       || (a[1] == b[1] && a[2] > 0.0)) ? 1 : 0;
}

int qd_gt ( const double *a, const double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] > b[1])
                       || (a[1] == b[1] && a[2] > b[2])
                       || (a[2] == b[2] && a[3] > b[3])) ? 1 : 0;
}

int qd_geq_qd_d ( const double *a, double b )
{
   return (a[0] > b || (a[0] == b && a[1] >= 0.0)) ? 1 : 0;
}

int qd_geq_qd_dd ( const double *a, const double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] > b[1])
                       || (a[1] == b[1] && a[2] >= 0.0)) ? 1 : 0;
}

int qd_geq ( const double *a, const double *b )
{
   return (a[0] > b[0] || (a[0] == b[0] && a[1] > b[1])
                       || (a[1] == b[1] && a[2] > b[2])
                       || (a[2] == b[2] && a[3] >= b[3])) ? 1 : 0;
}

/**************** Part II : operations from qd_real.cpp ****************/

/************************ divisions *********************************/

void qd_div_qd_d ( const double *a, double b, double *c )
{
/* Strategy: compute approximate quotient using high order doubles,
   and then correct it 3 times using the remainder (like long division).
*/
   double t0, t1, q0, q1, q2, q3;
   double dd_t[2];

   q0 = a[0]/b;                      /* approximate quotient */
   t0 = dd_two_prod(q0,b,&t1);       /* compute the remainder a - q0*b */
   dd_t[0] = t0; dd_t[1] = t1;       /* make double double for sub */
   qd_sub_qd_dd(a,dd_t,c);           /* c = a - dd_real(t0, t1); */
   q1 = c[0]/b;                      /* compute the first correction */
   t0 = dd_two_prod(q1,b,&t1);
   dd_t[0] = t0; dd_t[1] = t1;
   qd_sub_qd_dd(c,dd_t,c);           /* c -= dd_real(t0, t1); */
   q2 = c[0]/b;                      /* second correction to the quotient */
   t0 = dd_two_prod(q2,b,&t1);
   dd_t[0] = t0; dd_t[1] = t1;
   qd_sub_qd_dd(c,dd_t,c);           /* c -= dd_real(t0, t1); */
   q3 = c[0]/b;                      /* final correction to the quotient */
   qd_renorm4(&q0,&q1,&q2,&q3);
   c[0] = q0; c[1] = q1; c[2] = q2; c[3] = q3;
}

void qd_div_qd_dd ( const double *a, const double *b, double *c )
{
/* note: this is the accurate division */

   double q0, q1, q2, q3, q4;
   double qd_b[4];
   double acc[4];
   double r[4];

   qd_b[0] = b[0]; qd_b[1] = b[1]; qd_b[2] = 0.0; qd_b[3] = 0.0;
   q0 = a[0]/b[0];
   qd_mul_d(qd_b,q0,acc);      /* r = a - q0 * qd_b; */
   qd_sub(a,acc,r);
   q1 = r[0]/b[0];
   qd_mul_d(qd_b,q1,acc);
   qd_sub(r,acc,r);            /* r -= (q1 * qd_b); */
   q2 = r[0]/b[0];
   qd_mul_d(qd_b,q2,acc);
   qd_sub(r,acc,r);            /* r -= (q2 * qd_b); */
   q3 = r[0]/b[0];
   qd_mul_d(qd_b,q3,acc);
   qd_sub(r,acc,r);            /* r -= (q3 * qd_b); */
   q4 = r[0]/b[0];
   qd_renorm5(&q0,&q1,&q2,&q3,&q4);
   c[0] = q0; c[1] = q1; c[2] = q2; c[3] = q3;
}

void qd_div ( const double *a, const double *b, double *c )
{
/* note: this is the accurate division */

   double q0, q1, q2, q3, q4;
   double acc[4];
   double r[4];

   q0 = a[0]/b[0];
   qd_mul_d(b,q0,acc);
   qd_sub(a,acc,r);       /* r = a - (b * q0); */
   q1 = r[0]/b[0];
   qd_mul_d(b,q1,acc);
   qd_sub(r,acc,r);       /* r -= (b * q1); */
   q2 = r[0]/b[0];
   qd_mul_d(b,q2,acc);
   qd_sub(r,acc,r);       /* r -= (b * q2); */
   q3 = r[0]/b[0];
   qd_mul_d(b,q3,acc);
   qd_sub(r,acc,r);       /* r -= (b * q3); */
   q4 = r[0]/b[0];
   qd_renorm5(&q0,&q1,&q2,&q3,&q4);
   c[0] = q0; c[1] = q1; c[2] = q2; c[3] = q3;
}

void qd_div_d_qd ( double a, const double *b, double *c )
{
   double qd_a[4];

   qd_a[0] = a;   qd_a[1] = 0.0; 
   qd_a[2] = 0.0; qd_a[3] = 0.0;
   qd_div(qd_a,b,c);
}

void qd_div_dd_qd ( const double *a, const double *b, double *c )
{
   double qd_a[4];

   qd_a[0] = a[0]; qd_a[1] = a[1]; 
   qd_a[2] = 0.0;  qd_a[3] = 0.0;
   qd_div(qd_a,b,c);
}

/********************** squaring and power ***************************/

void qd_sqr ( const double *a, double *b )
{
/* quad-double ^ 2  = (x0 + x1 + x2 + x3) ^ 2
   = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2) + (2 x0 * x3 + 2 x1 * x2)
*/
   double p0, p1, p2, p3, p4, p5;
   double q0, q1, q2, q3;
   double s0, s1, t0, t1;

   p0 = dd_two_sqr(a[0],&q0);
   p1 = dd_two_prod(2.0*a[0],a[1],&q1);
   p2 = dd_two_prod(2.0*a[0],a[2],&q2);
   p3 = dd_two_sqr(a[1],&q3);
   p1 = dd_two_sum(q0,p1,&q0);
   q0 = dd_two_sum(q0,q1,&q1);
   p2 = dd_two_sum(p2,p3,&p3);
   s0 = dd_two_sum(q0,p2,&t0);
   s1 = dd_two_sum(q1,p3,&t1);
   s1 = dd_two_sum(s1,t0,&t0);
   t0 += t1;
   s1 = dd_quick_two_sum(s1,t0,&t0);
   p2 = dd_quick_two_sum(s0,s1,&t1);
   p3 = dd_quick_two_sum(t1,t0,&q0);
   p4 = 2.0*a[0]*a[3];
   p5 = 2.0*a[1]*a[2];
   p4 = dd_two_sum(p4,p5,&p5);
   q2 = dd_two_sum(q2,q3,&q3);
   t0 = dd_two_sum(p4,q2,&t1);
   t1 = t1 + p5 + q3;
   p3 = dd_two_sum(p3,t0,&p4);
   p4 = p4 + q0 + t1;
   qd_renorm5(&p0,&p1,&p2,&p3,&p4);
   b[0] = p0; b[1] = p1; b[2] = p2; b[3] = p3;
}

void qd_npwr ( const double *a, int n, double *b )
{
   if(n == 0)
   {
      b[0] = 1.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0;
   }
   else
   {
      int N = (n < 0) ? -n : n;   /* N = abs(n) */
      double r[4];                /* odd case multiplier */

      qd_copy(a,r);
      b[0] = 1.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0;
      if(N > 1)
      {                          /* use binary exponentiation */
         while(N > 0)
         {
            if (N % 2 == 1)      /* if odd multiply by r */
            {                    /* eventually N = 1, so this executes */
               qd_mul(b,r,b);
            }
            N /= 2;
            if(N > 0) qd_sqr(r,r); /* r = sqr(r); */
         }
      }
      else
      {
         qd_copy(r,b);
      }
      if(n < 0) qd_div_d_qd(1.0,b,b);
   }
}

void qd_ldexp ( double *a, int n )
{
   a[0] = ldexp(a[0],n); a[1] = ldexp(a[1],n);
   a[2] = ldexp(a[2],n); a[3] = ldexp(a[3],n);
}

/*********************** exp and log ******************************/

void qd_exp ( const double *a, double *b )
{
/* Strategy:  We first reduce the size of x by noting that

          exp(kr + m * log(2)) = 2^m * exp(r)^k

   where m and k are integers.  By choosing m appropriately
   we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
   evaluated using the familiar Taylor series.  Reducing the
   argument substantially speeds up the convergence.
*/
   const double k = ldexp(1.0,16);
   const double inv_k = 1.0/k;
   const double qd_e_0 =  2.718281828459045091e+00; /* exp(1)[0] */
   const double qd_e_1 =  1.445646891729250158e-16; /* exp(1)[1] */
   const double qd_e_2 = -2.127717108038176765e-33; /* exp(1)[2] */
   const double qd_e_3 =  1.515630159841218954e-49; /* exp(1)[3] */
   const double qd_log2_0 =  6.931471805599452862e-01; /* log(2)[0] */
   const double qd_log2_1 =  2.319046813846299558e-17; /* log(2)[1] */
   const double qd_log2_2 =  5.707708438416212066e-34; /* log(2)[2] */
   const double qd_log2_3 = -3.582432210601811423e-50; /* log(2)[3] */
   const double qd_eps = 1.21543267145725e-63;      /* 2^-209 */
   double m,tol;
   double r[4];
   double s[4];
   double p[4];
   double t[4];
   int i;
   double *inv_fac;
   double i_fac[60]; /* inverse factorials for Taylor expansion */
   i_fac[0] =  1.66666666666666657e-01;  i_fac[1] =  9.25185853854297066e-18;
   i_fac[2] =  5.13581318503262866e-34;  i_fac[3] =  2.85094902409834186e-50;
   i_fac[4] =  4.16666666666666644e-02;  i_fac[5] =  2.31296463463574266e-18;
   i_fac[6] =  1.28395329625815716e-34;  i_fac[7] =  7.12737256024585466e-51;
   i_fac[8] =  8.33333333333333322e-03;  i_fac[9] =  1.15648231731787138e-19;
   i_fac[10] =  1.60494162032269652e-36; i_fac[11] =  2.22730392507682967e-53;
   i_fac[12] =  1.38888888888888894e-03; i_fac[13] = -5.30054395437357706e-20;
   i_fac[14] = -1.73868675534958776e-36; i_fac[15] = -1.63335621172300840e-52;
   i_fac[16] =  1.98412698412698413e-04; i_fac[17] =  1.72095582934207053e-22;
   i_fac[18] =  1.49269123913941271e-40; i_fac[19] =  1.29470326746002471e-58;
   i_fac[20] =  2.48015873015873016e-05; i_fac[21] =  2.15119478667758816e-23;
   i_fac[22] =  1.86586404892426588e-41; i_fac[23] =  1.61837908432503088e-59;
   i_fac[24] =  2.75573192239858925e-06; i_fac[25] = -1.85839327404647208e-22;
   i_fac[26] =  8.49175460488199287e-39; i_fac[27] = -5.72661640789429621e-55;
   i_fac[28] =  2.75573192239858883e-07; i_fac[29] =  2.37677146222502973e-23;
   i_fac[30] = -3.26318890334088294e-40; i_fac[31] =  1.61435111860404415e-56;
   i_fac[32] =  2.50521083854417202e-08; i_fac[33] = -1.44881407093591197e-24;
   i_fac[34] =  2.04267351467144546e-41; i_fac[35] = -8.49632672007163175e-58;
   i_fac[36] =  2.08767569878681002e-09; i_fac[37] = -1.20734505911325997e-25;
   i_fac[38] =  1.70222792889287100e-42; i_fac[39] =  1.41609532150396700e-58;
   i_fac[40] =  1.60590438368216133e-10; i_fac[41] =  1.25852945887520981e-26;
   i_fac[42] = -5.31334602762985031e-43; i_fac[43] =  3.54021472597605528e-59;
   i_fac[44] =  1.14707455977297245e-11; i_fac[45] =  2.06555127528307454e-28;
   i_fac[46] =  6.88907923246664603e-45; i_fac[47] =  5.72920002655109095e-61;
   i_fac[48] =  7.64716373181981641e-13; i_fac[49] =  7.03872877733453001e-30;
   i_fac[50] = -7.82753927716258345e-48; i_fac[51] =  1.92138649443790242e-64;
   i_fac[52] =  4.77947733238738525e-14; i_fac[53] =  4.39920548583408126e-31;
   i_fac[54] = -4.89221204822661465e-49; i_fac[55] =  1.20086655902368901e-65;
   i_fac[56] =  2.81145725434552060e-15; i_fac[57] =  1.65088427308614326e-31;
   i_fac[58] = -2.87777179307447918e-50; i_fac[59] =  4.27110689256293549e-67;

   if(a[0] <= -709.0) 
   {
      b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; return;
   }
   if(a[0] >= 709.0) 
   {
      b[0] = -1.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; return;
   }
   if(qd_is_zero(a) == 1) 
   {
      b[0] = 1.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; return;
   }
   if(qd_is_one(a) == 1)
   {
      b[0] = qd_e_0; b[1] = qd_e_1; b[2] = qd_e_2; b[3] = qd_e_3; return;
   }
   m = floor(a[0] / qd_log2_0 + 0.5);
   /* r = mul_pwr2(a - qd_real::_log2 * m, inv_k); */
   r[0] = qd_log2_0; r[1] = qd_log2_1;
   r[2] = qd_log2_2; r[3] = qd_log2_3;
   qd_mul_d(r,m,r);
   qd_sub(a,r,r);
   qd_mul_pwr2(r,inv_k,r);
   qd_sqr(r,p);
   qd_mul_pwr2(p,0.5,s);   /* s = r + mul_pwr2(p,0.5); */
   qd_add(s,r,s);
   i = 0;
   tol = inv_k*qd_eps;
   do
   {
      qd_mul(r,p,p);              /* p *= r */
      inv_fac = &i_fac[4*(i++)];
      qd_mul(p,inv_fac,t);        /* t = p * inv_fac */
      qd_add(s,t,s);              /* s += t */
   } while(fabs(t[0]) > tol && i < 9);

   for(i=0; i<16; i++)     /* 16 times s = mul_pwr2(s,2.0) + sqr(s); */
   {
      qd_mul_pwr2(s,2.0,p);
      qd_sqr(s,t);
      qd_add(p,t,s);
   }
   qd_add_d(s,1.0,b);      /* s += 1.0; but then b stores the result */
   i = (int) m;
   qd_ldexp(b,i);
}

void qd_log ( const double *a, double *b )
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

     Three iterations are needed, since Newton's iteration
     approximately doubles the number of digits per iteration.
*/
   double acc[4];
   int i;

   if(qd_is_one(a) == 1)
   {
      b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; return;
   }
   if (a[0] <= 0.0)
   {
      printf("qd_log: argument is not positive");
      b[0] = -1.0; b[0] = 0.0; b[2] = 0.0; b[3] = 0.0; return;
   }
   b[0] = log(a[0]);              /* initial approximation */
   b[1] = 0.0; b[2] = 0.0; b[3] = 0.0;
   for(i=0; i<3; i++)
   {                              /* b = b + a*exp(-b) - 1 */
      qd_copy(b,acc);
      qd_minus(acc);              /* acc = -b */
      qd_exp(acc,acc);            /* acc = exp(-b) */
      qd_mul(a,acc,acc);          /* acc = a*exp(-b) */
      qd_add(b,acc,b);            /* b = b + a*exp(-b) */
      qd_sub_qd_d(b,1.0,b);       /* b = b + a*exp(-b) - 1 */
   }
}

void qd_log10 ( const double *a, double *b )
{
   const double qd_log10_0 =  2.302585092994045901e+00;
   const double qd_log10_1 = -2.170756223382249351e-16;
   const double qd_log10_2 = -9.984262454465776570e-33;
   const double qd_log10_3 = -4.023357454450206379e-49;
   double logten[4];
   
   logten[0] = qd_log10_0; logten[1] = qd_log10_1;
   logten[2] = qd_log10_2; logten[3] = qd_log10_3;
   qd_log(a,b);
   qd_div(b,logten,b);     /* b = log(a)/log(10) */
}

/********** Part III : input/output operations *****************/

int qd_read ( const char *s, double *a )
{
   const char *p = s;
   char ch;
   int sign = 0;
   int point = -1;         /* location of decimal point */
   int nd = 0;             /* number of digits read */
   int e = 0;              /* exponent */
   int done = 0;
   int nread;

   while(*p == ' ') p++;   /* skip any leading spaces */
   qd_real(0.0,a);
   while((done == 0) && (ch = *p) != '\0')
   {
      if(ch >= '0' && ch <= '9')    /* it is a digit */
      {
         int d = ch - '0';
         double cd = (double) d;

         qd_mul_d(a,10.0,a);        /* r *= 10.0; */
         qd_add_d(a,cd,a);          /* r += (double) d; */
         nd++;
      }
      else                          /* not a digit */
      {
         switch(ch)
         {
            case '.':
               if(point >= 0) return -1; /* already encountered point */
               point = nd; break;
            case '-':
            case '+':
               if(sign != 0 || nd > 0) return -1; /* already had sign */
               sign = (ch == '-') ? -1 : 1; break;
            case 'E':
            case 'e':
               nread = sscanf(p+1, "%d", &e);
               done = 1;
               if(nread != 1) return -1; /* reading of exponent failed */
               break;
            case ' ':
               done = 1; break;
            default:
               return -1;
         }
      }
      p++;
   }
   if(point >= 0) e -= (nd - point); /* adjust exponent for decimal point */
   if (e != 0)                       /* multiply the exponent */
   {
      double acc[4];

      qd_real(10.0,acc);
      qd_npwr(acc,e,acc);
      qd_mul(a,acc,a);               /* r *= (qd_real(10.0) ^ e); */
   }
   if(sign < 0) qd_minus(a);         /* qd = (sign < 0) ? -r : r; */

   return 0;
}

void qd_to_digits ( const double *a, char *s, int *expn, int precision )
{
   int D = precision + 1;  /* number of digits to compute */
   double acc[4];          /* for the r in the original code */
   double tmp[4];          /* to hold temporary results */
   int i,d,e;              /* e is exponent */

   qd_abs(a,acc);

   if(a[0] == 0.0)         /* a equals zero */
   {
      expn = 0;
      for(i = 0; i < precision; i++) s[i] = '0'; return;
   }
   e = (int) (floor(log10(fabs(a[0]))));  /* approximate exponent */
   if (e < -300)
   {
      qd_real(10.0,tmp);
      qd_npwr(tmp,300,tmp);
      qd_mul(acc,tmp,acc);             /* r *= dd_real(10.0) ^ 300; */
      qd_real(10.0,tmp);
      qd_npwr(tmp,e+300,tmp);
      qd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ (e + 300); */
   }
   else if(e > 300)
   {
      qd_ldexp(acc,-53);               /* r = ldexp(r, -53); */
      qd_real(10.0,tmp);
      qd_npwr(tmp,e,tmp);
      qd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ e; */
      qd_ldexp(acc,53);                /* r = ldexp(r, 53); */
   }
   else
   {
      qd_real(10.0,tmp);
      qd_npwr(tmp,e,tmp);
      qd_div(acc,tmp,acc);             /* r /= dd_real(10.0) ^ e; */
   }
   if(qd_gt_qd_d(acc,10.0) == 1)  /* fix exponent if we are off by one */
   {
      qd_div_qd_d(acc,10.0,acc);       /* r /= 10.0; */
      e++;
   }
   else if(qd_lt_qd_d(acc,1.0) == 1)   /* r < 1.0) */
   {
      qd_mul_d(acc,10.0,acc);          /* r *= 10.0; */
      e--;
   }
   if((qd_geq_qd_d(acc,10.0) == 1) || (qd_lt_qd_d(acc,1.0) == 1))
   {                 /* (r >= 10.0 || r < 1.0) */
      printf("qd_to_digits: cannot compute exponent"); return;
   }
   for(i=0; i<D; i++)                  /* extract the digits */
   {
      d = ((int) acc[0]);
      qd_sub_qd_d(acc,(double) d,acc); /* r -= d; */
      qd_mul_d(acc,10.0,acc);          /* r *= 10.0; */
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
      printf("qd_to_digits: nonpositive leading digit"); return;
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

void qd_to_string ( const double *a, char *s, int precision, int width,
                    int fixed, int showpos, int uppercase, char fill,
                    int *endstring )
{
/* note: nan and inf are ignored ...*/

   int sgn = 1;
   int i, e = 0;
   char *p = s;
   int cnt = 0; /* counts characters written to string */
   int off,d;
   double acc[4];

   if(qd_is_negative(a) == 1)
   {
      *(p++) = '-'; cnt++;
   }
   else if(showpos == 1)
   {
      *(p++) = '+'; cnt++;
   }
   else
      sgn = 0;

   if(qd_is_zero(a) == 1)
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
         qd_abs(a,acc);
         qd_log(acc,acc);
         qd_floor(acc,acc);
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
         char t[d+1];
         int j;

         qd_to_digits(a,t,&e,d);
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

void qd_write ( const double *a, int precision )
{
   char s[precision+10];
   int s_end;

   qd_to_string(a,s,precision,0,0,0,1,' ',&s_end);
   printf("%s",s);
}
