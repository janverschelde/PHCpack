/* file: triple_double.c */

/* This file contains the corresponding C code for the functions
   with prototypes declared in the triple_double.h file. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "double_double.h"
#include "triple_double.h"

/************************* normalizations ************************/

void td_fast_renorm
 ( double x0, double x1, double x2, double x3,
   double *r0, double *r1, double *r2 )
{
   double f0,f1,f2,f3,pr;
   int ptr;

   pr = dd_quick_two_sum(x2,x3,&f3);
   pr = dd_quick_two_sum(x1,pr,&f2);
   f0 = dd_quick_two_sum(x0,pr,&f1);
   if(f1 == 0.0)
   {
      pr = f0;
      ptr = 0;
      *r0 = dd_quick_two_sum(pr,f2,&pr);
   }
   else
   {
      *r0 = f0;
      pr = f1;
      ptr = 1;
      *r1 = dd_quick_two_sum(pr,f2,&pr);
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
      *r0 = dd_quick_two_sum(pr,f3,&pr);
   else if(ptr == 1)
      *r1 = dd_quick_two_sum(pr,f3,&pr);
   else
      *r2 = dd_quick_two_sum(pr,f3,&pr);

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

/****************************** copy *****************************/

void td_copy ( const double *a, double *b )
{
   b[0] = a[0];
   b[1] = a[1];
   b[2] = a[2];
}

/******************* addition and subtraction *********************/

void td_add ( const double *a, const double *b, double *c )
{
   double f0,f1,f2,f3,e;

   f2 = dd_two_sum(a[2],b[2],&f3);
   f1 = dd_two_sum(a[1],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;
   f0 = dd_two_sum(a[0],b[0],&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;

   td_fast_renorm(f0,f1,f2,f3,&c[0],&c[1],&c[2]);
}

void td_minus ( double *a )
{
   a[0] = -a[0];
   a[1] = -a[1];
   a[2] = -a[2];
}

void td_sub ( const double *a, const double *b, double *c )
{
   double minb[3];

   td_copy(b,minb);
   td_minus(minb);
   td_add(a,minb,c);
}

/**************  multiplication and division ***********************/

void td_mul ( const double *a, const double *b, double *c )
{
   double f0,f1,f2,f3,p,e;

   f3 = a[1]*b[2] + a[2]*b[1];
   f2 = dd_two_prod(a[0],b[2],&e);
   f3 += e;
   p = dd_two_prod(a[1],b[1],&e);
   f3 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 += e;
   p = dd_two_prod(a[2],b[0],&e);
   f3 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 += e;
   f1 = dd_two_prod(a[0],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;
   p = dd_two_prod(a[1],b[0],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;
   f1 = dd_two_sum(f1,p,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;
   f0 = dd_two_prod(a[0],b[0],&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;

   td_fast_renorm(f0,f1,f2,f3,&c[0],&c[1],&c[2]);
}

void td_mul_td_d ( const double *a, double b, double *c )
{
   double f0,f1,f2,f3,e;

   f3 = 0.0;
   f2 = dd_two_prod(a[2],b,&e);
   f3 += e;
   f1 = dd_two_prod(a[1],b,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;
   f0 = dd_two_prod(a[0],b,&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 += e;

   td_fast_renorm(f0,f1,f2,f3,&c[0],&c[1],&c[2]);
}

void td_div ( const double *a, const double *b, double *c )
{
   double acc[3];
   double q0,q1,q2,q3;

   q0 = a[0]/b[0];
   td_mul_td_d(b,q0,acc);
   td_sub(a,acc,c);

   q1 = c[0]/b[0];
   td_mul_td_d(b,q1,acc);
   td_sub(c,acc,c);

   q2 = c[0]/b[0];
   td_mul_td_d(b,q2,acc);
   td_sub(c,acc,c);

   q3 = c[0]/b[0];

   td_fast_renorm(q0,q1,q2,q3,&c[0],&c[1],&c[2]);
}

/**************************** random *****************************/

void td_random ( double *x )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double dd_eps = 4.930380657631324e-32; // 2^(-104)
   const double first = ((double) rand())/RAND_MAX;
   const double second = ((double) rand())/RAND_MAX;
   const double third = ((double) rand())/RAND_MAX;

   td_fast_renorm(first,eps*second,dd_eps*third,0.0,&x[0],&x[1],&x[2]);
}

/************************ basic output *********************************/

void td_write_doubles ( const double* x )
{
   printf("  hi = %21.14e\n",x[0]);
   printf("  mi = %21.14e\n",x[1]);
   printf("  lo = %21.14e\n",x[2]);
}
