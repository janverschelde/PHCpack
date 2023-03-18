/* file: deca_double.c */

/* This file contains the corresponding C code for the functions
   with prototypes declared in the deca_double.h file. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "double_double.h"
#include "deca_double.h"

/************************* normalizations ************************/

void da_renorm10
 ( double f0, double f1, double f2, double f3, double f4, double f5,
   double f6, double f7, double f8, double f9, double f10, double *pr,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9 )
{
   int ptr;

   if(f1 == 0.0)
   {
      *pr = f0;
      ptr = 0;
      *r0 = dd_quick_two_sum(*pr,f2,pr);
   }
   else
   {
      *r0 = f0;
      *pr = f1;
      ptr = 1;
      *r1 = dd_quick_two_sum(*pr,f2,pr);
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
      *r0 = dd_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f3,pr);
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
      *r0 = dd_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f4,pr);
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
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f5,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f6,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
      else if(ptr == 5)
      {
         *pr = *r5;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f7,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
      else if(ptr == 5)
      {
         *pr = *r5;
      }
      else if(ptr == 6)
      {
         *pr = *r6;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f8,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
      else if(ptr == 5)
      {
         *pr = *r5;
      }
      else if(ptr == 6)
      {
         *pr = *r6;
      }
      else if(ptr == 7)
      {
         *pr = *r7;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f9,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
      else if(ptr == 5)
      {
         *pr = *r5;
      }
      else if(ptr == 6)
      {
         *pr = *r6;
      }
      else if(ptr == 7)
      {
         *pr = *r7;
      }
      else if(ptr == 8)
      {
         *pr = *r8;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)
   {
      *r0 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f10,pr);
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
      else if(ptr == 4)
      {
         *pr = *r4;
      }
      else if(ptr == 5)
      {
         *pr = *r5;
      }
      else if(ptr == 6)
      {
         *pr = *r6;
      }
      else if(ptr == 7)
      {
         *pr = *r7;
      }
      else if(ptr == 8)
      {
         *pr = *r8;
      }
      else if(ptr == 9)
      {
         *pr = *r9;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if((ptr < 10) && (*pr != 0.0))
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
      else if(ptr == 4)
      {
         *r4 = *pr;
      }
      else if(ptr == 5)
      {
         *r5 = *pr;
      }
      else if(ptr == 6)
      {
         *r6 = *pr;
      }
      else if(ptr == 7)
      {
         *r7 = *pr;
      }
      else if(ptr == 8)
      {
         *r8 = *pr;
      }
      else if(ptr == 9)
      {
         *r9 = *pr;
      }
      ptr = ptr + 1;
   }
   if(ptr < 1)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0; *r0 = 0.0;
   }
   else if(ptr < 2)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0;
   }
   else if(ptr < 3)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0;
   }
   else if(ptr < 4)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0;
   }
   else if(ptr < 5)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
   }
   else if(ptr < 6)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0;
   }
   else if(ptr < 7)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0;
   }
   else if(ptr < 8)
   {
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0;
   }
   else if(ptr < 9)
   {
      *r9 = 0.0; *r8 = 0.0;
   }
   else if(ptr < 10)
   {
      *r9 = 0.0;
   }
}

void da_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr;

   pr = dd_quick_two_sum(x9,x10,&f10);
   pr = dd_quick_two_sum(x8,pr,&f9);
   pr = dd_quick_two_sum(x7,pr,&f8);
   pr = dd_quick_two_sum(x6,pr,&f7);
   pr = dd_quick_two_sum(x5,pr,&f6);
   pr = dd_quick_two_sum(x4,pr,&f5);
   pr = dd_quick_two_sum(x3,pr,&f4);
   pr = dd_quick_two_sum(x2,pr,&f3);
   pr = dd_quick_two_sum(x1,pr,&f2);
   f0 = dd_quick_two_sum(x0,pr,&f1);

   da_renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,&pr,
               r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
}

void da_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr;

   pr = dd_two_sum(x9,y,&f10);
   pr = dd_two_sum(x8,pr,&f9);
   pr = dd_two_sum(x7,pr,&f8);
   pr = dd_two_sum(x6,pr,&f7);
   pr = dd_two_sum(x5,pr,&f6);
   pr = dd_two_sum(x4,pr,&f5);
   pr = dd_two_sum(x3,pr,&f4);
   pr = dd_two_sum(x2,pr,&f3);
   pr = dd_two_sum(x1,pr,&f2);
   f0 = dd_two_sum(x0,pr,&f1);

   da_renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,&pr,
               r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
}

/****************************** copy *****************************/

void da_copy ( const double *a, double *b )
{
   b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; b[3] = a[3]; b[4] = a[4];
   b[5] = a[5]; b[6] = a[6]; b[7] = a[7]; b[8] = a[8]; b[9] = a[9];
}

/******************* addition and subtraction *********************/

void da_add ( const double *a, const double *b, double *c )
{
   // ALGORITHM : baileyAdd_fast<10,10,10> generated by CAMPARY.
   
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e;

   f10 = 0.0;
   f9 = dd_two_sum(a[9],b[9],&e);
   f10 += e;
   f8 = dd_two_sum(a[8],b[8],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(a[7],b[7],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(a[6],b[6],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(a[5],b[5],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_sum(a[4],b[4],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_sum(a[3],b[3],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f2 = dd_two_sum(a[2],b[2],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f1 = dd_two_sum(a[1],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f0 = dd_two_sum(a[0],b[0],&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;

   da_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                  &c[0],&c[1],&c[2],&c[3],&c[4],
                  &c[5],&c[6],&c[7],&c[8],&c[9]);
}

void da_add_da_d ( const double *a, double b, double *c )
{
   da_renorm_add1
      (a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],b,
       &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9]);
}

void da_minus ( double *a )
{
   a[0] = -a[0]; a[1] = -a[1]; a[2] = -a[2]; a[3] = -a[3]; a[4] = -a[4];
   a[5] = -a[5]; a[6] = -a[6]; a[7] = -a[7]; a[8] = -a[8]; a[9] = -a[9];
}

void da_sub ( const double *a, const double *b, double *c )
{
   double minb[10];

   da_copy(b,minb);
   da_minus(minb);
   da_add(a,minb,c);
}

/**************  multiplications and division ***********************/

void da_mul ( const double *a, const double *b, double *c )
{
   // ALGORITHM :baileyMul_fast<10,10,10> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,p,e;

   f10 =  a[1]*b[9];
   f10 += a[2]*b[8];
   f10 += a[3]*b[7];
   f10 += a[4]*b[6];
   f10 += a[5]*b[5];
   f10 += a[6]*b[4];
   f10 += a[7]*b[3];
   f10 += a[8]*b[2];
   f10 += a[9]*b[1];
   f9 = dd_two_prod(a[0],b[9],&e);
   f10 += e;
   p = dd_two_prod(a[1],b[8],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[7],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[6],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[5],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[5],b[4],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[6],b[3],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[7],b[2],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[8],b[1],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   p = dd_two_prod(a[9],b[0],&e);
   f10 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 += e;
   f8 = dd_two_prod(a[0],b[8],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[7],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[6],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[5],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[4],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[5],b[3],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[6],b[2],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[7],b[1],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[8],b[0],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_prod(a[0],b[7],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[6],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[5],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[4],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[3],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[5],b[2],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[6],b[1],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[7],b[0],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_prod(a[0],b[6],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[5],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[4],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[3],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[2],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[5],b[1],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[6],b[0],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_prod(a[0],b[5],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[4],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[3],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[2],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[1],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[5],b[0],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_prod(a[0],b[4],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[3],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[2],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[1],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[4],b[0],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_prod(a[0],b[3],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[2],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[1],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[3],b[0],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f2 = dd_two_prod(a[0],b[2],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[1],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[2],b[0],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f1 = dd_two_prod(a[0],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   p = dd_two_prod(a[1],b[0],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f1 = dd_two_sum(f1,p,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f0 = dd_two_prod(a[0],b[0],&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;

   da_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                  &c[0],&c[1],&c[2],&c[3],&c[4],
                  &c[5],&c[6],&c[7],&c[8],&c[9]);
}

void da_mul_da_d ( const double *a, double b, double *c )
{
   // ALGORITHM : baileyMul_fast<10,1,10> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e;

   f10 = 0.0;
   f9 = dd_two_prod(a[9],b,&e);
   f10 += e;
   f8 = dd_two_prod(a[8],b,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f7 = dd_two_prod(a[7],b,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f6 = dd_two_prod(a[6],b,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f5 = dd_two_prod(a[5],b,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f4 = dd_two_prod(a[4],b,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f3 = dd_two_prod(a[3],b,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f2 = dd_two_prod(a[2],b,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f1 = dd_two_prod(a[1],b,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;
   f0 = dd_two_prod(a[0],b,&e);
   f1 = dd_two_sum(f1,e,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 += e;

   da_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                  &c[0],&c[1],&c[2],&c[3],&c[4],
                  &c[5],&c[6],&c[7],&c[8],&c[9]);
}

void da_div ( const double *a, const double *b, double *c )
{
   double acc[10];
   double q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10;

   q0 = a[0]/b[0]; da_mul_da_d(b,q0,acc); da_sub(a,acc,c);
   q1 = c[0]/b[0]; da_mul_da_d(b,q1,acc); da_sub(c,acc,c);
   q2 = c[0]/b[0]; da_mul_da_d(b,q2,acc); da_sub(c,acc,c);
   q3 = c[0]/b[0]; da_mul_da_d(b,q3,acc); da_sub(c,acc,c);
   q4 = c[0]/b[0]; da_mul_da_d(b,q4,acc); da_sub(c,acc,c);
   q5 = c[0]/b[0]; da_mul_da_d(b,q5,acc); da_sub(c,acc,c);
   q6 = c[0]/b[0]; da_mul_da_d(b,q6,acc); da_sub(c,acc,c);
   q7 = c[0]/b[0]; da_mul_da_d(b,q7,acc); da_sub(c,acc,c);
   q8 = c[0]/b[0]; da_mul_da_d(b,q8,acc); da_sub(c,acc,c);
   q9 = c[0]/b[0]; da_mul_da_d(b,q9,acc); da_sub(c,acc,c);
   q10 = c[0]/b[0];

   da_fast_renorm
      (q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,
       &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9]);
}

/******************* random number generator ***************************/

void da_random ( double *x )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double eps4 = eps*eps3;
   const double eps5 = eps*eps4;
   const double eps6 = eps*eps5;
   const double eps7 = eps*eps6;
   const double eps8 = eps*eps7;
   const double eps9 = eps*eps8;
   const double r0 = ((double) rand())/RAND_MAX;
   const double r1 = ((double) rand())/RAND_MAX;
   const double r2 = ((double) rand())/RAND_MAX;
   const double r3 = ((double) rand())/RAND_MAX;
   const double r4 = ((double) rand())/RAND_MAX;
   const double r5 = ((double) rand())/RAND_MAX;
   const double r6 = ((double) rand())/RAND_MAX;
   const double r7 = ((double) rand())/RAND_MAX;
   const double r8 = ((double) rand())/RAND_MAX;
   const double r9 = ((double) rand())/RAND_MAX;

   x[0] = r0; x[1] = 0.0; x[2] = 0.0; x[3] = 0.0; x[4] = 0.0;
   x[5] = 0.0; x[6] = 0.0; x[7] = 0.0; x[8] = 0.0; x[9] = 0.0;
   da_add_da_d(x,r1*eps,x);
   da_add_da_d(x,r2*eps2,x);
   da_add_da_d(x,r3*eps3,x);
   da_add_da_d(x,r4*eps4,x);
   da_add_da_d(x,r5*eps5,x);
   da_add_da_d(x,r6*eps6,x);
   da_add_da_d(x,r7*eps7,x);
   da_add_da_d(x,r8*eps8,x);
   da_add_da_d(x,r9*eps9,x);
}

/************************ basic output *********************************/

void da_write_doubles ( const double *x )
{
   printf("  R thumb  = %21.14e",x[0]);
   printf("  R index  = %21.14e\n",x[1]);
   printf("  R middle = %21.14e",x[2]);
   printf("  R ring   = %21.14e\n",x[3]);
   printf("  R pink   = %21.14e",x[4]);
   printf("  L thumb  = %21.14e\n",x[5]);
   printf("  L index  = %21.14e",x[6]);
   printf("  L middle = %21.14e\n",x[7]);
   printf("  L ring   = %21.14e",x[8]);
   printf("  L pink   = %21.14e\n",x[9]);
}
