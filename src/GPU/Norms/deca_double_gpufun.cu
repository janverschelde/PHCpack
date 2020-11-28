// The file deca_double_gpufun.cpp defines the code for the functions
// specified in deca_double_gpufun.h

#include "double_double_gpufun.h"
#include "octo_double_gpufun.h"
#include "deca_double_gpufun.h"

/************************* renormalizations ************************/

__device__ __forceinline__ void dag_renorm10
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
      *r0 = ddg_quick_two_sum(*pr,f2,pr);
   }
   else
   {
      *r0 = f0;
      *pr = f1;
      ptr = 1;
      *r1 = ddg_quick_two_sum(*pr,f2,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f3,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f4,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f5,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f6,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f7,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f8,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f9,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f9,pr);
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
      *r0 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f10,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f10,pr);
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

__device__ __forceinline__ void dag_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr;

   pr = ddg_quick_two_sum(x9,x10,&f10);
   pr = ddg_quick_two_sum(x8,pr,&f9);
   pr = ddg_quick_two_sum(x7,pr,&f8);
   pr = ddg_quick_two_sum(x6,pr,&f7);
   pr = ddg_quick_two_sum(x5,pr,&f6);
   pr = ddg_quick_two_sum(x4,pr,&f5);
   pr = ddg_quick_two_sum(x3,pr,&f4);
   pr = ddg_quick_two_sum(x2,pr,&f3);
   pr = ddg_quick_two_sum(x1,pr,&f2);
   f0 = ddg_quick_two_sum(x0,pr,&f1);

   dag_renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,&pr,
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
}

__device__ __forceinline__ void dag_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr;

   pr = ddg_two_sum(x9,y,&f10);
   pr = ddg_two_sum(x8,pr,&f9);
   pr = ddg_two_sum(x7,pr,&f8);
   pr = ddg_two_sum(x6,pr,&f7);
   pr = ddg_two_sum(x5,pr,&f6);
   pr = ddg_two_sum(x4,pr,&f5);
   pr = ddg_two_sum(x3,pr,&f4);
   pr = ddg_two_sum(x2,pr,&f3);
   pr = ddg_two_sum(x1,pr,&f2);
   f0 = ddg_two_sum(x0,pr,&f1);

   dag_renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,&pr,
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
}

/************************ copy and abs *******************************/

__device__ __forceinline__ void dag_copy
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double *b_rtb, double *b_rix, double *b_rmi, double *b_rrg, double *b_rpk,
   double *b_ltb, double *b_lix, double *b_lmi, double *b_lrg, double *b_lpk )
{
   *b_rtb = a_rtb;
   *b_rix = a_rix;
   *b_rmi = a_rmi;
   *b_rrg = a_rrg;
   *b_rpk = a_rpk;
   *b_ltb = a_ltb;
   *b_lix = a_lix;
   *b_lmi = a_lmi;
   *b_lrg = a_lrg;
   *b_lpk = a_lpk;
}

__device__ __forceinline__ void dag_abs
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double *b_rtb, double *b_rix, double *b_rmi, double *b_rrg, double *b_rpk,
   double *b_ltb, double *b_lix, double *b_lmi, double *b_lrg, double *b_lpk )
{
   if(a_rtb < 0.0)
   {
      *b_rtb = -a_rtb;
      *b_rix = -a_rix;
      *b_rmi = -a_rmi;
      *b_rrg = -a_rrg;
      *b_rpk = -a_rpk;
      *b_ltb = -a_ltb;
      *b_lix = -a_lix;
      *b_lmi = -a_lmi;
      *b_lrg = -a_lrg;
      *b_lpk = -a_lpk;
   }
   else
   {
      *b_rtb = a_rtb;
      *b_rix = a_rix;
      *b_rmi = a_rmi;
      *b_rrg = a_rrg;
      *b_rpk = a_rpk;
      *b_ltb = a_ltb;
      *b_lix = a_lix;
      *b_lmi = a_lmi;
      *b_lrg = a_lrg;
      *b_lpk = a_lpk;
   }
}

/****************** additions and substractions ************************/

__device__ __forceinline__ void dag_add
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b_rtb, double b_rix, double b_rmi, double b_rrg, double b_rpk,
   double b_ltb, double b_lix, double b_lmi, double b_lrg, double b_lpk,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   // ALGORITHM : baileyAddg_fast<10,10,10> generated by CAMPARY.
   
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e;

   f10 = 0.0;
   f9 = ddg_two_sum(a_lpk,b_lpk,&e);
   f10 += e;
   f8 = ddg_two_sum(a_lrg,b_lrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(a_lmi,b_lmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(a_lix,b_lix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(a_ltb,b_ltb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(a_rpk,b_rpk,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(a_rrg,b_rrg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(a_rmi,b_rmi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_sum(a_rix,b_rix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f0 = ddg_two_sum(a_rtb,b_rtb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;

   dag_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                   c_rtb,c_rix,c_rmi,c_rrg,c_rpk,
                   c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
}

__device__ __forceinline__ void dag_inc
 ( double *a_rtb, double *a_rix, double *a_rmi, double *a_rrg, double *a_rpk,
   double *a_ltb, double *a_lix, double *a_lmi, double *a_lrg, double *a_lpk,
   double b_rtb, double b_rix, double b_rmi, double b_rrg, double b_rpk,
   double b_ltb, double b_lix, double b_lmi, double b_lrg, double b_lpk )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e;

   f10 = 0.0;
   f9 = ddg_two_sum(*a_lpk,b_lpk,&e);
   f10 += e;
   f8 = ddg_two_sum(*a_lrg,b_lrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(*a_lmi,b_lmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(*a_lix,b_lix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(*a_ltb,b_ltb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(*a_rpk,b_rpk,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(*a_rrg,b_rrg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(*a_rmi,b_rmi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_sum(*a_rix,b_rix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f0 = ddg_two_sum(*a_rtb,b_rtb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;

   dag_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                   a_rtb,a_rix,a_rmi,a_rrg,a_rpk,
                   a_ltb,a_lix,a_lmi,a_lrg,a_lpk);
}

__device__ __forceinline__ void dag_inc_d
 ( double *a_rtb, double *a_rix, double *a_rmi, double *a_rrg, double *a_rpk,
   double *a_ltb, double *a_lix, double *a_lmi, double *a_lrg, double *a_lpk,
   double b )
{
   dag_renorm_add1
      (*a_rtb,*a_rix,*a_rmi,*a_rrg,*a_rpk,*a_ltb,*a_lix,*a_lmi,*a_lrg,*a_lpk,
       b,a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,a_lrg,a_lpk);
}

__device__ __forceinline__ void dag_minus
 ( double *a_rtb, double *a_rix, double *a_rmi, double *a_rrg, double *a_rpk,
   double *a_ltb, double *a_lix, double *a_lmi, double *a_lrg, double *a_lpk )
{
   *a_rtb = -(*a_rtb);
   *a_rix = -(*a_rix);
   *a_rmi = -(*a_rmi);
   *a_rrg = -(*a_rrg);
   *a_rpk = -(*a_rpk);
   *a_ltb = -(*a_ltb);
   *a_lix = -(*a_lix);
   *a_lmi = -(*a_lmi);
   *a_lrg = -(*a_lrg);
   *a_lpk = -(*a_lpk);
}

__device__ __forceinline__ void dag_sub
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b_rtb, double b_rix, double b_rmi, double b_rrg, double b_rpk,
   double b_ltb, double b_lix, double b_lmi, double b_lrg, double b_lpk,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   dag_copy(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,b_ltb,b_lix,b_lmi,b_lrg,b_lpk,
            c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
   dag_minus(c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
   dag_inc(c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk,
           a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,a_lrg,a_lpk);
}

/***************** multiplications and division ********************/

__device__ __forceinline__ void dag_mul_pwr2
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   *c_rtb = a_rtb*b;
   *c_rix = a_rix*b;
   *c_rmi = a_rmi*b;
   *c_rrg = a_rrg*b;
   *c_rpk = a_rpk*b;
   *c_ltb = a_ltb*b;
   *c_lix = a_lix*b;
   *c_lmi = a_lmi*b;
   *c_lrg = a_lrg*b;
   *c_lpk = a_lpk*b;
}

__device__ __forceinline__ void dag_mul
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b_rtb, double b_rix, double b_rmi, double b_rrg, double b_rpk,
   double b_ltb, double b_lix, double b_lmi, double b_lrg, double b_lpk,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   // ALGORITHM :baileyMul_fast<10,10,10> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,p,e;

   f10 =  a_rix*b_lpk;
   f10 += a_rmi*b_lrg;
   f10 += a_rrg*b_lmi;
   f10 += a_rpk*b_lix;
   f10 += a_ltb*b_ltb;
   f10 += a_lix*b_rpk;
   f10 += a_lmi*b_rrg;
   f10 += a_lrg*b_rmi;
   f10 += a_lpk*b_rix;
   f9 = ddg_two_prod(a_rtb,b_lpk,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_lrg,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_lmi,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_lix,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_ltb,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,b_rpk,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,b_rrg,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,b_rmi,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lrg,b_rix,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lpk,b_rtb,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   f8 = ddg_two_prod(a_rtb,b_lrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_lmi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_lix,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_ltb,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_rpk,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,b_rrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,b_rmi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,b_rix,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lrg,b_rtb,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_prod(a_rtb,b_lmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_lix,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_ltb,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_rpk,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_rrg,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,b_rmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,b_rix,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,b_rtb,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_prod(a_rtb,b_lix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_ltb,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_rpk,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_rrg,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_rmi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,b_rix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,b_rtb,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_prod(a_rtb,b_ltb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_rpk,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_rrg,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_rmi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_rix,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,b_rtb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_prod(a_rtb,b_rpk,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_rrg,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_rmi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_rix,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,b_rtb,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_prod(a_rtb,b_rrg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_rmi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_rix,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,b_rtb,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_prod(a_rtb,b_rmi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_rix,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,b_rtb,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_prod(a_rtb,b_rix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,b_rtb,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f0 = ddg_two_prod(a_rtb,b_rtb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;

   dag_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                   c_rtb,c_rix,c_rmi,c_rrg,c_rpk,
                   c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
}

__device__ __forceinline__ void dag_sqr
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,p,e;

   f10 =  a_rix*a_lpk;
   f10 += a_rmi*a_lrg;
   f10 += a_rrg*a_lmi;
   f10 += a_rpk*a_lix;
   f10 += a_ltb*a_ltb;
   f10 += a_lix*a_rpk;
   f10 += a_lmi*a_rrg;
   f10 += a_lrg*a_rmi;
   f10 += a_lpk*a_rix;
   f9 = ddg_two_prod(a_rtb,a_lpk,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_lrg,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_lmi,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_lix,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_ltb,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,a_rpk,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,a_rrg,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,a_rmi,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lrg,a_rix,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   p = ddg_two_prod(a_lpk,a_rtb,&e);
   f10 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 += e;
   f8 = ddg_two_prod(a_rtb,a_lrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_lmi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_lix,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_ltb,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_rpk,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,a_rrg,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,a_rmi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,a_rix,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lrg,a_rtb,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_prod(a_rtb,a_lmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_lix,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_ltb,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_rpk,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_rrg,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,a_rmi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,a_rix,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lmi,a_rtb,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_prod(a_rtb,a_lix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_ltb,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_rpk,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_rrg,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_rmi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,a_rix,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_lix,a_rtb,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_prod(a_rtb,a_ltb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_rpk,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_rrg,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_rmi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_rix,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_ltb,a_rtb,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_prod(a_rtb,a_rpk,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_rrg,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_rmi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_rix,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rpk,a_rtb,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_prod(a_rtb,a_rrg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_rmi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_rix,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rrg,a_rtb,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_prod(a_rtb,a_rmi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_rix,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rmi,a_rtb,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_prod(a_rtb,a_rix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   p = ddg_two_prod(a_rix,a_rtb,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f0 = ddg_two_prod(a_rtb,a_rtb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;

   dag_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                   c_rtb,c_rix,c_rmi,c_rrg,c_rpk,
                   c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
}

__device__ __forceinline__ void dag_mul_da_d
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   // ALGORITHM : baileyMul_fast<10,1,10> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e;

   f10 = 0.0;
   f9 = ddg_two_prod(a_lpk,b,&e);
   f10 += e;
   f8 = ddg_two_prod(a_lrg,b,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f7 = ddg_two_prod(a_lmi,b,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f6 = ddg_two_prod(a_lix,b,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f5 = ddg_two_prod(a_ltb,b,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f4 = ddg_two_prod(a_rpk,b,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f3 = ddg_two_prod(a_rrg,b,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f2 = ddg_two_prod(a_rmi,b,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f1 = ddg_two_prod(a_rix,b,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;
   f0 = ddg_two_prod(a_rtb,b,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 += e;

   dag_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                   c_rtb,c_rix,c_rmi,c_rrg,c_rpk,
                   c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
}

__device__ __forceinline__ void dag_div
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double b_rtb, double b_rix, double b_rmi, double b_rrg, double b_rpk,
   double b_ltb, double b_lix, double b_lmi, double b_lrg, double b_lpk,
   double *c_rtb, double *c_rix, double *c_rmi, double *c_rrg, double *c_rpk,
   double *c_ltb, double *c_lix, double *c_lmi, double *c_lrg, double *c_lpk )
{
   double acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk;
   double acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk;
   double q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10;

   q0 = a_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q0,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,a_lrg,a_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q1 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q1,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q2 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q2,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q3 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q3,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q4 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q4,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q5 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q5,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q6 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q6,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q7 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q7,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q8 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q8,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q9 = *c_rtb/b_rtb;
   dag_mul_da_d(b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk,q9,
                &acc_rtb,&acc_rix,&acc_rmi,&acc_rrg,&acc_rpk,
                &acc_ltb,&acc_lix,&acc_lmi,&acc_lrg,&acc_lpk);
   dag_sub(*c_rtb,*c_rix,*c_rmi,*c_rrg,*c_rpk,
           *c_ltb,*c_lix,*c_lmi,*c_lrg,*c_lpk,
           acc_rtb,acc_rix,acc_rmi,acc_rrg,acc_rpk,
           acc_ltb,acc_lix,acc_lmi,acc_lrg,acc_lpk,
           c_rtb,c_rix,c_rmi,c_rrg,c_rpk,c_ltb,c_lix,c_lmi,c_lrg,c_lpk);

   q10 = *c_rtb/b_rtb;

   dag_fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,
                   c_rtb,c_rix,c_rmi,c_rrg,c_rpk,
                   c_ltb,c_lix,c_lmi,c_lrg,c_lpk);
}

/***************************** square root *****************************/

__device__ __forceinline__ void dag_sqrt
 ( double a_rtb, double a_rix, double a_rmi, double a_rrg, double a_rpk,
   double a_ltb, double a_lix, double a_lmi, double a_lrg, double a_lpk,
   double *b_rtb, double *b_rix, double *b_rmi, double *b_rrg, double *b_rpk,
   double *b_ltb, double *b_lix, double *b_lmi, double *b_lrg, double *b_lpk )
{
   double z_rtb,z_rix,z_rmi,z_rrg,z_rpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk;

   odg_sqrt(a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,
            b_rtb,b_rix,b_rmi,b_rrg,b_rpk,b_ltb,b_lix,b_lmi);

   *b_lrg = 0.0; *b_lpk = 0.0;

   dag_sqr(*b_rtb,*b_rix,*b_rmi,*b_rrg,*b_rpk,
           *b_ltb,*b_lix,*b_lmi,*b_lrg,*b_lpk,
           &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
           &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
   dag_inc(&z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
           &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk,
           a_rtb,a_rix,a_rmi,a_rrg,a_rpk,
           a_ltb,a_lix,a_lmi,a_lrg,a_lpk);
   dag_div(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk,
           *b_rtb,*b_rix,*b_rmi,*b_rrg,*b_rpk,
           *b_ltb,*b_lix,*b_lmi,*b_lrg,*b_lpk,
           &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
           &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
   dag_mul_pwr2(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,
                z_ltb,z_lix,z_lmi,z_lrg,z_lpk,0.5,
                b_rtb,b_rix,b_rmi,b_rrg,b_rpk,
                b_ltb,b_lix,b_lmi,b_lrg,b_lpk);
}
