// The file penta_double_gpufun.cu defines the code for the functions
// specified in penta_double_gpufun.h

#include "double_double_gpufun.h"
#include "penta_double_gpufun.h"

/************************** renormalizations **************************/

__device__ __forceinline__ void pdg_renorm5
 ( double f0, double f1, double f2, double f3, double f4, double f5,
   double *pr, double *r0, double *r1, double *r2, double *r3, double *r4 )
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
   if((ptr < 5) && (*pr != 0.0))
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
      ptr = ptr + 1;
   }
   if(ptr < 1)
   {
      *r4 = 0.0; *r3 = 0.0; *r2 = 0.0; *r1 = 0.0; *r0 = 0.0;
   }
   else if(ptr < 2)
   {
      *r4 = 0.0; *r3 = 0.0; *r2 = 0.0; *r1 = 0.0;
   }
   else if(ptr < 3)
   {
      *r4 = 0.0; *r3 = 0.0; *r2 = 0.0;
   }
   else if(ptr < 4)
   {
      *r4 = 0.0; *r3 = 0.0;
   }
   else if(ptr < 5)
   {
      *r4 = 0.0;
   }
}

__device__ __forceinline__ void pdg_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double *r0, double *r1, double *r2, double *r3, double *r4 )
{
   double f0,f1,f2,f3,f4,f5,pr;

   pr = ddg_quick_two_sum(x4,x5,&f5);
   pr = ddg_quick_two_sum(x3,pr,&f4);
   pr = ddg_quick_two_sum(x2,pr,&f3);
   pr = ddg_quick_two_sum(x1,pr,&f2);
   f0 = ddg_quick_two_sum(x0,pr,&f1);

   pdg_renorm5(f0,f1,f2,f3,f4,f5,&pr,r0,r1,r2,r3,r4);
}

__device__ __forceinline__ void pdg_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4 )
{
   double f0,f1,f2,f3,f4,f5,pr;

   pr = ddg_two_sum(x4,y,&f5);
   pr = ddg_two_sum(x3,pr,&f4);
   pr = ddg_two_sum(x2,pr,&f3);
   pr = ddg_two_sum(x1,pr,&f2);
   f0 = ddg_two_sum(x0,pr,&f1);

   pdg_renorm5(f0,f1,f2,f3,f4,f5,&pr,r0,r1,r2,r3,r4);
}

/************************ copy and abs *******************************/

__device__ __forceinline__ void pdg_copy
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double *b_tb, double *b_ix, double *b_mi, double *b_rg, double *b_pk )
{
   *b_tb = a_tb;
   *b_ix = a_ix;
   *b_mi = a_mi;
   *b_rg = a_rg;
   *b_pk = a_pk;
}

__device__ __forceinline__ void pdg_abs
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double *b_tb, double *b_ix, double *b_mi, double *b_rg, double *b_pk )
{
   if(a_tb < 0.0)
   {
      *b_tb = -a_tb;
      *b_ix = -a_ix;
      *b_mi = -a_mi;
      *b_rg = -a_rg;
      *b_pk = -a_pk;
   }
   else
   {
      *b_tb = a_tb;
      *b_ix = a_ix;
      *b_mi = a_mi;
      *b_rg = a_rg;
      *b_pk = a_pk;
   }
}

/****************** additions and substractions ************************/

__device__ __forceinline__ void pdg_add
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b_tb, double b_ix, double b_mi, double b_rg, double b_pk,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk)
{
   // ALGORITHM : baileyAdd_fast<5,5,5> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,e;

   f5 = 0.0;
   f4 = ddg_two_sum(a_pk,b_pk,&e);
   f5 += e;
   f3 = ddg_two_sum(a_rg,b_rg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(a_mi,b_mi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_sum(a_ix,b_ix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f0 = ddg_two_sum(a_tb,b_tb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;

   pdg_fast_renorm(f0,f1,f2,f3,f4,f5,c_tb,c_ix,c_mi,c_rg,c_pk);
}

__device__ __forceinline__ void pdg_inc
 ( double *a_tb, double *a_ix, double *a_mi, double *a_rg, double *a_pk,
   double b_tb, double b_ix, double b_mi, double b_rg, double b_pk )
{
   double f0,f1,f2,f3,f4,f5,e;

   f5 = 0.0;
   f4 = ddg_two_sum(*a_pk,b_pk,&e);
   f5 += e;
   f3 = ddg_two_sum(*a_rg,b_rg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(*a_mi,b_mi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_sum(*a_ix,b_ix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f0 = ddg_two_sum(*a_tb,b_tb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;

   pdg_fast_renorm(f0,f1,f2,f3,f4,f5,a_tb,a_ix,a_mi,a_rg,a_pk);
}

__device__ __forceinline__ void pdg_inc_d
 ( double *a_tb, double *a_ix, double *a_mi, double *a_rg, double *a_pk,
   double b )
{
   pdg_renorm_add1(*a_tb,*a_ix,*a_mi,*a_rg,*a_pk,b,a_tb,a_ix,a_mi,a_rg,a_pk);
}

__device__ __forceinline__ void pdg_minus
 ( double *a_tb, double *a_ix, double *a_mi, double *a_rg, double *a_pk )
{
   *a_tb = -(*a_tb);
   *a_ix = -(*a_ix);
   *a_mi = -(*a_mi);
   *a_rg = -(*a_rg);
   *a_pk = -(*a_pk);
}

__device__ __forceinline__ void pdg_sub
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b_tb, double b_ix, double b_mi, double b_rg, double b_pk,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   pdg_copy(b_tb,b_ix,b_mi,b_rg,b_pk,c_tb,c_ix,c_mi,c_rg,c_pk);
   pdg_minus(c_tb,c_ix,c_mi,c_rg,c_pk);
   pdg_inc(c_tb,c_ix,c_mi,c_rg,c_pk,a_tb,a_ix,a_mi,a_rg,a_pk);
}

/***************** multiplications and division ********************/

__device__ __forceinline__ void pdg_mul_pwr2
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   *c_tb = a_tb*b;
   *c_ix = a_ix*b;
   *c_mi = a_mi*b;
   *c_rg = a_rg*b;
   *c_pk = a_pk*b;
}

__device__ __forceinline__ void pdg_mul
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b_tb, double b_ix, double b_mi, double b_rg, double b_pk,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   // ALGORITHM : baileyMul_fast<5,5,5> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,p,e;

   f5 =  a_ix*b_pk;
   f5 += a_mi*b_rg;
   f5 += a_rg*b_mi;
   f5 += a_pk*b_ix;

   f4 = ddg_two_prod(a_tb,b_pk,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,b_rg,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,b_mi,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_rg,b_ix,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_pk,b_tb,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   f3 = ddg_two_prod(a_tb,b_rg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,b_mi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,b_ix,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_rg,b_tb,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_prod(a_tb,b_mi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,b_ix,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,b_tb,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_prod(a_tb,b_ix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,b_tb,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f0 = ddg_two_prod(a_tb,b_tb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;

   pdg_fast_renorm(f0,f1,f2,f3,f4,f5,c_tb,c_ix,c_mi,c_rg,c_pk);
}

__device__ __forceinline__ void pdg_sqr
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   double f0,f1,f2,f3,f4,f5,p,e;

   f5 =  a_ix*a_pk;
   f5 += a_mi*a_rg;
   f5 += a_rg*a_mi;
   f5 += a_pk*a_ix;

   f4 = ddg_two_prod(a_tb,a_pk,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,a_rg,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,a_mi,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_rg,a_ix,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   p = ddg_two_prod(a_pk,a_tb,&e);
   f5 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 += e;
   f3 = ddg_two_prod(a_tb,a_rg,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,a_mi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,a_ix,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_rg,a_tb,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_prod(a_tb,a_mi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,a_ix,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_mi,a_tb,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_prod(a_tb,a_ix,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   p = ddg_two_prod(a_ix,a_tb,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f0 = ddg_two_prod(a_tb,a_tb,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;

   pdg_fast_renorm(f0,f1,f2,f3,f4,f5,c_tb,c_ix,c_mi,c_rg,c_pk);
}

__device__ __forceinline__ void pdg_mul_pd_d
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   // ALGORITHM : baileyMul_fast<5,1,5>

   double f0,f1,f2,f3,f4,f5,e;

   f5 = 0.0;
   f4 = ddg_two_prod(a_pk,b,&e);
   f5 += e;
   f3 = ddg_two_prod(a_rg,b,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f2 = ddg_two_prod(a_mi,b,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f1 = ddg_two_prod(a_ix,b,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;
   f0 = ddg_two_prod(a_tb,b,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 += e;

   pdg_fast_renorm(f0,f1,f2,f3,f4,f5,c_tb,c_ix,c_mi,c_rg,c_pk);
}

__device__ __forceinline__ void pdg_div
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double b_tb, double b_ix, double b_mi, double b_rg, double b_pk,
   double *c_tb, double *c_ix, double *c_mi, double *c_rg, double *c_pk )
{
   double acc_tb,acc_ix,acc_mi,acc_rg,acc_pk;
   double q0,q1,q2,q3,q4,q5;

   q0 = a_tb/b_tb;
   pdg_mul_pd_d(b_tb,b_ix,b_mi,b_rg,b_pk,q0,
                &acc_tb,&acc_ix,&acc_mi,&acc_rg,&acc_pk);
   pdg_sub(a_tb,a_ix,a_mi,a_rg,a_pk,acc_tb,acc_ix,acc_mi,acc_rg,acc_pk,
           c_tb,c_ix,c_mi,c_rg,c_pk);

   q1 = *c_tb/b_tb;
   pdg_mul_pd_d(b_tb,b_ix,b_mi,b_rg,b_pk,q1,
                &acc_tb,&acc_ix,&acc_mi,&acc_rg,&acc_pk);
   pdg_sub(*c_tb,*c_ix,*c_mi,*c_rg,*c_pk,acc_tb,acc_ix,acc_mi,acc_rg,acc_pk,
           c_tb,c_ix,c_mi,c_rg,c_pk);

   q2 = *c_tb/b_tb;
   pdg_mul_pd_d(b_tb,b_ix,b_mi,b_rg,b_pk,q2,
                &acc_tb,&acc_ix,&acc_mi,&acc_rg,&acc_pk);
   pdg_sub(*c_tb,*c_ix,*c_mi,*c_rg,*c_pk,acc_tb,acc_ix,acc_mi,acc_rg,acc_pk,
           c_tb,c_ix,c_mi,c_rg,c_pk);

   q3 = *c_tb/b_tb;
   pdg_mul_pd_d(b_tb,b_ix,b_mi,b_rg,b_pk,q3,
                &acc_tb,&acc_ix,&acc_mi,&acc_rg,&acc_pk);
   pdg_sub(*c_tb,*c_ix,*c_mi,*c_rg,*c_pk,acc_tb,acc_ix,acc_mi,acc_rg,acc_pk,
           c_tb,c_ix,c_mi,c_rg,c_pk);

   q4 = *c_tb/b_tb;
   pdg_mul_pd_d(b_tb,b_ix,b_mi,b_rg,b_pk,q4,
                &acc_tb,&acc_ix,&acc_mi,&acc_rg,&acc_pk);
   pdg_sub(*c_tb,*c_ix,*c_mi,*c_rg,*c_pk,acc_tb,acc_ix,acc_mi,acc_rg,acc_pk,
           c_tb,c_ix,c_mi,c_rg,c_pk);

   q5 = *c_tb/b_tb;

   pdg_fast_renorm(q0,q1,q2,q3,q4,q5,c_tb,c_ix,c_mi,c_rg,c_pk);
}

/***************************** square root *****************************/

__device__ __forceinline__ void pdg_sqrt
 ( double a_tb, double a_ix, double a_mi, double a_rg, double a_pk,
   double *b_tb, double *b_ix, double *b_mi, double *b_rg, double *b_pk )
{
   double z_tb,z_ix,z_mi,z_rg,z_pk;

   ddg_sqrt(a_tb,a_ix,b_tb,b_ix);
   pdg_sqr(*b_tb,*b_ix,0.0,0.0,0.0,&z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
   pdg_inc(&z_tb,&z_ix,&z_mi,&z_rg,&z_pk,a_tb,a_ix,a_mi,a_rg,a_pk);
   pdg_div(z_tb,z_ix,z_mi,z_rg,z_pk,*b_tb,*b_ix,0.0,0.0,0.0,
           &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
   pdg_mul_pwr2(z_tb,z_ix,z_mi,z_rg,z_pk,0.5,b_tb,b_ix,b_mi,b_rg,b_pk);
   pdg_sqr(*b_tb,*b_ix,*b_mi,*b_rg,*b_pk,&z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
   pdg_inc(&z_tb,&z_ix,&z_mi,&z_rg,&z_pk,a_tb,a_ix,a_mi,a_rg,a_pk);
   pdg_div(z_tb,z_ix,z_mi,z_rg,z_pk,*b_tb,*b_ix,*b_mi,*b_rg,*b_pk,
           &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
   pdg_mul_pwr2(z_tb,z_ix,z_mi,z_rg,z_pk,0.5,b_tb,b_ix,b_mi,b_rg,b_pk);
}
