// The file hexa_double_gpufun.cu defines the code for the functions
// specified in hexa_double_gpufun.h

#include "double_double_gpufun.h"
#include "quad_double_gpufun.h"
#include "octo_double_gpufun.h"
#include "hexa_double_gpufun.h"

/************************* renormalizations ************************/

__device__ __forceinline__ void hdg_renorm16
 ( double f0, double f1, double f2, double f3, double f4, double f5,
   double f6, double f7, double f8, double f9, double f10, double f11,
   double f12, double f13, double f14, double f15, double f16, double *pr, 
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
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
   if(*pr == 0.0)    // first test on pr
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
   if(ptr == 0)     // first test on ptr
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
   if(*pr == 0.0)   // second test on pr
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
   if(ptr == 0)    // second test on ptr
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
   if(*pr == 0.0)  // third test on pr
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
   if(ptr == 0)   // third test on ptr
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
   if(*pr == 0.0)   // fourth test on pr
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
   if(ptr == 0)   // fourth test on ptr
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
   if(*pr == 0.0)   // fifth test on pr
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
   if(ptr == 0)  // fifth test on ptr
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
   if(*pr == 0.0)  // sixth test on pr
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
   if(ptr == 0)  // sixth test on ptr
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
   if(*pr == 0.0)  // seventh test on pr
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
   if(ptr == 0)  // seventh test on ptr
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
   if(*pr == 0.0)  // eighth test on pr
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
   if(ptr == 0)  // eighth test on ptr
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
   if(*pr == 0.0)  // nineth test on pr
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
   if(ptr == 0)  // nineth test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f11,pr);
   }
   if(*pr == 0.0)  // tenth test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)  // tenth test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 11)
   {
      *r11 = ddg_quick_two_sum(*pr,f12,pr);
   }
   if(*pr == 0.0)  // eleventh test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
      else if(ptr == 11)
      {
         *pr = *r11;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)  // eleventh test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 11)
   {
      *r11 = ddg_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 12)
   {
      *r12 = ddg_quick_two_sum(*pr,f13,pr);
   }
   if(*pr == 0.0)  // twelveth test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
      else if(ptr == 11)
      {
         *pr = *r11;
      }
      else if(ptr == 12)
      {
         *pr = *r12;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)  // twelveth test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 11)
   {
      *r11 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 12)
   {
      *r12 = ddg_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 13)
   {
      *r13 = ddg_quick_two_sum(*pr,f14,pr);
   }
   if(*pr == 0.0)  // thirteenth test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
      else if(ptr == 11)
      {
         *pr = *r11;
      }
      else if(ptr == 12)
      {
         *pr = *r12;
      }
      else if(ptr == 13)
      {
         *pr = *r13;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)  // thirteenth test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 11)
   {
      *r11 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 12)
   {
      *r12 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 13)
   {
      *r13 = ddg_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 14)
   {
      *r14 = ddg_quick_two_sum(*pr,f15,pr);
   }
   if(*pr == 0.0)  // forteenth test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
      else if(ptr == 11)
      {
         *pr = *r11;
      }
      else if(ptr == 12)
      {
         *pr = *r12;
      }
      else if(ptr == 13)
      {
         *pr = *r13;
      }
      else if(ptr == 14)
      {
         *pr = *r14;
      }
   }
   else
   {
      ptr = ptr + 1;
   }
   if(ptr == 0)  // forteenth test on ptr
   {
      *r0 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 8)
   {
      *r8 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 9)
   {
      *r9 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 10)
   {
      *r10 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 11)
   {
      *r11 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 12)
   {
      *r12 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 13)
   {
      *r13 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 14)
   {
      *r14 = ddg_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 15)
   {
      *r15 = ddg_quick_two_sum(*pr,f16,pr);
   }
   if(*pr == 0.0)  // fifteenth test on pr
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
      else if(ptr == 10)
      {
         *pr = *r10;
      }
      else if(ptr == 11)
      {
         *pr = *r11;
      }
      else if(ptr == 12)
      {
         *pr = *r12;
      }
      else if(ptr == 13)
      {
         *pr = *r13;
      }
      else if(ptr == 14)
      {
         *pr = *r14;
      }
      else if(ptr == 15)
      {
         *pr = *r15;
      }
   }
   else
   {
      ptr = ptr + 1;
   }

// beginning of the end ...

   if((ptr < 16) && (*pr != 0.0))
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
      else if(ptr == 10)
      {
         *r10 = *pr;
      }
      else if(ptr == 11)
      {
         *r11 = *pr;
      }
      else if(ptr == 12)
      {
         *r12 = *pr;
      }
      else if(ptr == 13)
      {
         *r13 = *pr;
      }
      else if(ptr == 14)
      {
         *r14 = *pr;
      }
      else if(ptr == 15)
      {
         *r15 = *pr;
      }
      ptr = ptr + 1;
   }
   if(ptr < 1)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0; *r0 = 0.0;
   }
   else if(ptr < 2)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0; *r1 = 0.0;
   }
   else if(ptr < 3)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0; *r2 = 0.0;
   }
   else if(ptr < 4)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
      *r3 = 0.0;
   }
   else if(ptr < 5)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
   }
   else if(ptr < 6)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0; *r5 = 0.0;
   }
   else if(ptr < 7)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0; *r6 = 0.0;
   }
   else if(ptr < 8)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0; *r7 = 0.0;
   }
   else if(ptr < 9)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0; *r8 = 0.0;
   }
   else if(ptr < 10)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
      *r9 = 0.0;
   }
   else if(ptr < 11)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0; *r10 = 0.0;
   }
   else if(ptr < 12)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0; *r11 = 0.0;
   }
   else if(ptr < 13)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0; *r12 = 0.0;
   }
   else if(ptr < 14)
   {
      *r15 = 0.0; *r14 = 0.0; *r13 = 0.0;
   }
   else if(ptr < 15)
   {
      *r15 = 0.0; *r14 = 0.0;
   }
   else if(ptr < 16)
   {
      *r15 = 0.0;
   }
}

__device__ __forceinline__ void hdg_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double x16,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr;

   pr = ddg_quick_two_sum(x15,x16,&f16);
   pr = ddg_quick_two_sum(x14,pr,&f15);
   pr = ddg_quick_two_sum(x13,pr,&f14);
   pr = ddg_quick_two_sum(x12,pr,&f13);
   pr = ddg_quick_two_sum(x11,pr,&f12);
   pr = ddg_quick_two_sum(x10,pr,&f11);
   pr = ddg_quick_two_sum(x9,pr,&f10);
   pr = ddg_quick_two_sum(x8,pr,&f9);
   pr = ddg_quick_two_sum(x7,pr,&f8);
   pr = ddg_quick_two_sum(x6,pr,&f7);
   pr = ddg_quick_two_sum(x5,pr,&f6);
   pr = ddg_quick_two_sum(x4,pr,&f5);
   pr = ddg_quick_two_sum(x3,pr,&f4);
   pr = ddg_quick_two_sum(x2,pr,&f3);
   pr = ddg_quick_two_sum(x1,pr,&f2);
   f0 = ddg_quick_two_sum(x0,pr,&f1);

   hdg_renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,&pr,
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
}

__device__ __forceinline__ void hdg_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr;

   pr = ddg_two_sum(x15,y,&f16);
   pr = ddg_two_sum(x14,pr,&f15);
   pr = ddg_two_sum(x13,pr,&f14);
   pr = ddg_two_sum(x12,pr,&f13);
   pr = ddg_two_sum(x11,pr,&f12);
   pr = ddg_two_sum(x10,pr,&f11);
   pr = ddg_two_sum(x9,pr,&f10);
   pr = ddg_two_sum(x8,pr,&f9);
   pr = ddg_two_sum(x7,pr,&f8);
   pr = ddg_two_sum(x6,pr,&f7);
   pr = ddg_two_sum(x5,pr,&f6);
   pr = ddg_two_sum(x4,pr,&f5);
   pr = ddg_two_sum(x3,pr,&f4);
   pr = ddg_two_sum(x2,pr,&f3);
   pr = ddg_two_sum(x1,pr,&f2);
   f0 = ddg_two_sum(x0,pr,&f1);

   hdg_renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,&pr,
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
}

/*********************** copy and abs *****************************/

__device__ __forceinline__ void hdg_copy 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
   *b_hihihihi = a_hihihihi; *b_lohihihi = a_lohihihi;
   *b_hilohihi = a_hilohihi; *b_lolohihi = a_lolohihi;
   *b_hihilohi = a_hihilohi; *b_lohilohi = a_lohilohi;
   *b_hilolohi = a_hilolohi; *b_lololohi = a_lololohi;
   *b_hihihilo = a_hihihilo; *b_lohihilo = a_lohihilo;
   *b_hilohilo = a_hilohilo; *b_lolohilo = a_lolohilo;
   *b_hihilolo = a_hihilolo; *b_lohilolo = a_lohilolo;
   *b_hilololo = a_hilololo; *b_lolololo = a_lolololo;
}

__device__ __forceinline__ void hdg_abs 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
   if(a_hihihihi < 0.0)
   {
      *b_hihihihi = -a_hihihihi; *b_lohihihi = -a_lohihihi;
      *b_hilohihi = -a_hilohihi; *b_lolohihi = -a_lolohihi;
      *b_hihilohi = -a_hihilohi; *b_lohilohi = -a_lohilohi;
      *b_hilolohi = -a_hilolohi; *b_lololohi = -a_lololohi;
      *b_hihihilo = -a_hihihilo; *b_lohihilo = -a_lohilohi;
      *b_hilohilo = -a_hilohilo; *b_lolohilo = -a_lolohilo;
      *b_hihilolo = -a_hihilolo; *b_lohilolo = -a_lohilolo;
      *b_hilololo = -a_hilololo; *b_lolololo = -a_lolololo;
   }
   else
   {
      *b_hihihihi = a_hihihihi; *b_lohihihi = a_lohihihi;
      *b_hilohihi = a_hilohihi; *b_lolohihi = a_lolohihi;
      *b_hihilohi = a_hihilohi; *b_lohilohi = a_lohilohi;
      *b_hilolohi = a_hilolohi; *b_lololohi = a_lololohi;
      *b_hihihilo = a_hihihilo; *b_lohihilo = a_lohilohi;
      *b_hilohilo = a_hilohilo; *b_lolohilo = a_lolohilo;
      *b_hihilolo = a_hihilolo; *b_lohilolo = a_lohilolo;
      *b_hilololo = a_hilololo; *b_lolololo = a_lolololo;
   }
}

/******************* addition and subtraction *********************/

__device__ __forceinline__ void hdg_add 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   // ALGORITHM : baileyAdd_fast<16,16,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,e;

   f16 = 0.0;
   f15 = ddg_two_sum(a_lolololo,b_lolololo,&e);
   f16 += e;
   f14 = ddg_two_sum(a_hilololo,b_hilololo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(a_lohilolo,b_lohilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(a_hihilolo,b_hihilolo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(a_lolohilo,b_lolohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(a_hilohilo,b_hilohilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(a_lohihilo,b_lohihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(a_hihihilo,b_hihihilo,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(a_lololohi,b_lololohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(a_hilolohi,b_hilolohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(a_lohilohi,b_lohilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(a_hihilohi,b_hihilohi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(a_lolohihi,b_lolohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(a_hilohihi,b_hilohihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f1 = ddg_two_sum(a_lohihihi,b_lohihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f0 = ddg_two_sum(a_hihihihi,b_hihihihi,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   hdg_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                   c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
                   c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
                   c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
                   c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
}

__device__ __forceinline__ void hdg_inc
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,e;

   f16 = 0.0;
   f15 = ddg_two_sum(*a_lolololo,b_lolololo,&e);
   f16 += e;
   f14 = ddg_two_sum(*a_hilololo,b_hilololo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(*a_lohilolo,b_lohilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(*a_hihilolo,b_hihilolo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(*a_lolohilo,b_lolohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(*a_hilohilo,b_hilohilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(*a_lohihilo,b_lohihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(*a_hihihilo,b_hihihilo,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(*a_lololohi,b_lololohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(*a_hilolohi,b_hilolohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(*a_lohilohi,b_lohilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(*a_hihilohi,b_hihilohi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(*a_lolohihi,b_lolohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(*a_hilohihi,b_hilohihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f1 = ddg_two_sum(*a_lohihihi,b_lohihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f0 = ddg_two_sum(*a_hihihihi,b_hihihihi,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   hdg_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                   a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
                   a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
                   a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
                   a_hihilolo,a_lohilolo,a_hilololo,a_lolololo);
}

__device__ __forceinline__ void hdg_inc_d
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b )
{
   hdg_renorm_add1(*a_hihihihi,*a_lohihihi,*a_hilohihi,*a_lolohihi,
                   *a_hihilohi,*a_lohilohi,*a_hilolohi,*a_lololohi,
                   *a_hihihilo,*a_lohihilo,*a_hilohilo,*a_lolohilo,
                   *a_hihilolo,*a_lohilolo,*a_hilololo,*a_lolololo,b,
                   a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
                   a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
                   a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
                   a_hihilolo,a_lohilolo,a_hilololo,a_lolololo);
}

__device__ __forceinline__ void hdg_dec
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo )
{
   const double mbhihihihi = -b_hihihihi;
   const double mblohihihi = -b_lohihihi;
   const double mbhilohihi = -b_hilohihi;
   const double mblolohihi = -b_lolohihi;
   const double mbhihilohi = -b_hihilohi;
   const double mblohilohi = -b_lohilohi;
   const double mbhilolohi = -b_hilolohi;
   const double mblololohi = -b_lololohi;
   const double mbhihihilo = -b_hihihilo;
   const double mblohihilo = -b_lohihilo;
   const double mbhilohilo = -b_hilohilo;
   const double mblolohilo = -b_lolohilo;
   const double mbhihilolo = -b_hihilolo;
   const double mblohilolo = -b_lohilolo;
   const double mbhilololo = -b_hilololo;
   const double mblolololo = -b_lolololo;

   hdg_inc(a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
           a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
           a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
           a_hihilolo,a_lohilolo,a_hilololo,a_lolololo,
           mbhihihihi,mblohihihi,mbhilohihi,mblolohihi,
           mbhihilohi,mblohilohi,mbhilolohi,mblololohi,
           mbhihihilo,mblohihilo,mbhilohilo,mblolohilo,
           mbhihilolo,mblohilolo,mbhilololo,mblolololo);
}

__device__ __forceinline__ void hdg_minus 
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo )
{
   *a_hihihihi = -(*a_hihihihi); *a_lohihihi = -(*a_lohihihi);
   *a_hilohihi = -(*a_hilohihi); *a_lolohihi = -(*a_lolohihi);
   *a_hihilohi = -(*a_hihilohi); *a_lohilohi = -(*a_lohilohi);
   *a_hilolohi = -(*a_hilolohi); *a_lololohi = -(*a_lololohi);
   *a_hihihilo = -(*a_hihihilo); *a_lohihilo = -(*a_lohihilo);
   *a_hilohilo = -(*a_hilohilo); *a_lolohilo = -(*a_lolohilo);
   *a_hihilolo = -(*a_hihilolo); *a_lohilolo = -(*a_lohilolo);
   *a_hilololo = -(*a_hilololo); *a_lolololo = -(*a_lolololo);
}

__device__ __forceinline__ void hdg_sub 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   hdg_copy(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
            b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
            b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
            b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,
            c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
            c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
            c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
            c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
   hdg_minus(c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
             c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
             c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
             c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
   hdg_inc(c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo,
           a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
           a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
           a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
           a_hihilolo,a_lohilolo,a_hilololo,a_lolololo);
}

/**************  multiplications and division ***********************/

__device__ __forceinline__ void hdg_mul_pwr2
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   *c_hihihihi = a_hihihihi*b;
   *c_lohihihi = a_lohihihi*b;
   *c_hilohihi = a_hilohihi*b;
   *c_lolohihi = a_lolohihi*b;
   *c_hihilohi = a_hihilohi*b;
   *c_lohilohi = a_lohilohi*b;
   *c_hilolohi = a_hilolohi*b;
   *c_lololohi = a_lololohi*b;
   *c_hihihilo = a_hihihilo*b;
   *c_lohihilo = a_lohihilo*b;
   *c_hilohilo = a_hilohilo*b;
   *c_lolohilo = a_lolohilo*b;
   *c_hihilolo = a_hihilolo*b;
   *c_lohilolo = a_lohilolo*b;
   *c_hilololo = a_hilololo*b;
   *c_lolololo = a_lolololo*b;
}

__device__ __forceinline__ void hdg_mul
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   // ALGORITHM : baileyMul_fast<16,16,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,p,e;

   f16 =  a_lohihihi*b_lolololo;
   f16 += a_hilohihi*b_hilololo;
   f16 += a_lolohihi*b_lohilolo;
   f16 += a_hihilohi*b_hihilolo;
   f16 += a_lohilohi*b_lolohilo;
   f16 += a_hilolohi*b_hilohilo;
   f16 += a_lololohi*b_lohihilo;
   f16 += a_hihihilo*b_hihihilo;
   f16 += a_lohihilo*b_lololohi;
   f16 += a_hilohilo*b_hilolohi;
   f16 += a_lolohilo*b_lohilohi;
   f16 += a_hihilolo*b_hihilohi;
   f16 += a_lohilolo*b_lolohihi;
   f16 += a_hilololo*b_hilohihi;
   f16 += a_lolololo*b_lohihihi;

   f15 = ddg_two_prod(a_hihihihi,b_lolololo,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hilololo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lohilolo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hihilolo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lolohilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hilohilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_lohihilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_hihihilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_lololohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_hilolohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_lohilohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,b_hihilohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,b_lolohihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,b_hilohihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilololo,b_lohihihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolololo,b_hihihihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;

   f14 = ddg_two_prod(a_hihihihi,b_hilololo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lohilolo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hihilolo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lolohilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hilohilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_lohihilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_hihihilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_lololohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_hilolohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_lohilohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_hihilohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,b_lolohihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,b_hilohihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,b_lohihihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilololo,b_hihihihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f13 = ddg_two_prod(a_hihihihi,b_lohilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hihilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lolohilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hilohilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lohihilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hihihilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_lololohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_hilolohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_lohilohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_hihilohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_lolohihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,b_hilohihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,b_lohihihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,b_hihihihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f12 = ddg_two_prod(a_hihihihi,b_hihilolo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lolohilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hilohilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lohihilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hihihilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_lololohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_hilolohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_lohilohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_hihilohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_lolohihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_hilohihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,b_lohihihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,b_hihihihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f11 = ddg_two_prod(a_hihihihi,b_lolohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hilohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lohihilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hihihilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lololohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hilolohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_lohilohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_hihilohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_lolohihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_hilohihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_lohihihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,b_hihihihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f10 = ddg_two_prod(a_hihihihi,b_hilohilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lohihilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hihihilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lololohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hilolohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_lohilohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_hihilohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_lolohihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_hilohihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_lohihihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,b_hihihihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f9 = ddg_two_prod(a_hihihihi,b_lohihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hihihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lololohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hilolohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lohilohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hihilohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_lolohihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_hilohihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_lohihihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,b_hihihihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f8 = ddg_two_prod(a_hihihihi,b_hihihilo,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lololohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hilolohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lohilohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hihilohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_lolohihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_hilohihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_lohihihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,b_hihihihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f7 = ddg_two_prod(a_hihihihi,b_lololohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hilolohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lohilohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hihilohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lolohihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hilohihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_lohihihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,b_hihihihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f6 = ddg_two_prod(a_hihihihi,b_hilolohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lohilohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hihilohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lolohihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hilohihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_lohihihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,b_hihihihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f5 = ddg_two_prod(a_hihihihi,b_lohilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hihilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lolohihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hilohihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_lohihihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,b_hihihihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f4 = ddg_two_prod(a_hihihihi,b_hihilohi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lolohihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hilohihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_lohihihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,b_hihihihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f3 = ddg_two_prod(a_hihihihi,b_lolohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hilohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_lohihihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,b_hihihihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f2 = ddg_two_prod(a_hihihihi,b_hilohihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_lohihihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,b_hihihihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f1 = ddg_two_prod(a_hihihihi,b_lohihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,b_hihihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f0 = ddg_two_prod(a_hihihihi,b_hihihihi,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   hdg_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                   c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
                   c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
                   c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
                   c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
}

__device__ __forceinline__ void hdg_sqr
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,p,e;

   f16 =  a_lohihihi*a_lolololo;
   f16 += a_hilohihi*a_hilololo;
   f16 += a_lolohihi*a_lohilolo;
   f16 += a_hihilohi*a_hihilolo;
   f16 += a_lohilohi*a_lolohilo;
   f16 += a_hilolohi*a_hilohilo;
   f16 += a_lololohi*a_lohihilo;
   f16 += a_hihihilo*a_hihihilo;
   f16 += a_lohihilo*a_lololohi;
   f16 += a_hilohilo*a_hilolohi;
   f16 += a_lolohilo*a_lohilohi;
   f16 += a_hihilolo*a_hihilohi;
   f16 += a_lohilolo*a_lolohihi;
   f16 += a_hilololo*a_hilohihi;
   f16 += a_lolololo*a_lohihihi;

   f15 = ddg_two_prod(a_hihihihi,a_lolololo,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hilololo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lohilolo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hihilolo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lolohilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hilohilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_lohihilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_hihihilo,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_lololohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_hilolohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_lohilohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,a_hihilohi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,a_lolohihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,a_hilohihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_hilololo,a_lohihihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;
   p = ddg_two_prod(a_lolololo,a_hihihihi,&e);
   f16 += e;
   f15 = ddg_two_sum(f15,p,&e);
   f16 += e;

   f14 = ddg_two_prod(a_hihihihi,a_hilololo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lohilolo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hihilolo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lolohilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hilohilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_lohihilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_hihihilo,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_lololohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_hilolohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_lohilohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_hihilohi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,a_lolohihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,a_hilohihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,a_lohihihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilololo,a_hihihihi,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f14 = ddg_two_sum(f14,p,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f13 = ddg_two_prod(a_hihihihi,a_lohilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hihilolo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lolohilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hilohilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lohihilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hihihilo,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_lololohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_hilolohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_lohilohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_hihilohi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_lolohihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,a_hilohihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,a_lohihihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilolo,a_hihihihi,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_sum(f13,p,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f12 = ddg_two_prod(a_hihihihi,a_hihilolo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lolohilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hilohilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lohihilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hihihilo,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_lololohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_hilolohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_lohilohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_hihilohi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_lolohihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_hilohihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,a_lohihihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilolo,a_hihihihi,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_sum(f12,p,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f11 = ddg_two_prod(a_hihihihi,a_lolohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hilohilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lohihilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hihihilo,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lololohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hilolohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_lohilohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_hihilohi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_lolohihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_hilohihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_lohihihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohilo,a_hihihihi,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_sum(f11,p,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f10 = ddg_two_prod(a_hihihihi,a_hilohilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lohihilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hihihilo,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lololohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hilolohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_lohilohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_hihilohi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_lolohihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_hilohihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_lohihihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohilo,a_hihihihi,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_sum(f10,p,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f9 = ddg_two_prod(a_hihihihi,a_lohihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hihihilo,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lololohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hilolohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lohilohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hihilohi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_lolohihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_hilohihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_lohihihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihilo,a_hihihihi,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_sum(f9,p,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f8 = ddg_two_prod(a_hihihihi,a_hihihilo,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lololohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hilolohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lohilohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hihilohi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_lolohihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_hilohihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_lohihihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihihilo,a_hihihihi,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_sum(f8,p,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f7 = ddg_two_prod(a_hihihihi,a_lololohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hilolohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lohilohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hihilohi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lolohihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hilohihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_lohihihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lololohi,a_hihihihi,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_sum(f7,p,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f6 = ddg_two_prod(a_hihihihi,a_hilolohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lohilohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hihilohi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lolohihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hilohihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_lohihihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilolohi,a_hihihihi,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_sum(f6,p,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f5 = ddg_two_prod(a_hihihihi,a_lohilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hihilohi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lolohihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hilohihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_lohihihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohilohi,a_hihihihi,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_sum(f5,p,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f4 = ddg_two_prod(a_hihihihi,a_hihilohi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lolohihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hilohihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_lohihihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hihilohi,a_hihihihi,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_sum(f4,p,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f3 = ddg_two_prod(a_hihihihi,a_lolohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hilohihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_lohihihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lolohihi,a_hihihihi,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_sum(f3,p,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f2 = ddg_two_prod(a_hihihihi,a_hilohihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_lohihihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_hilohihi,a_hihihihi,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_sum(f2,p,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f1 = ddg_two_prod(a_hihihihi,a_lohihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   p = ddg_two_prod(a_lohihihi,a_hihihihi,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f1 = ddg_two_sum(f1,p,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   f0 = ddg_two_prod(a_hihihihi,a_hihihihi,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   hdg_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                   c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
                   c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
                   c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
                   c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
}

__device__ __forceinline__ void hdg_mul_hd_d
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   // ALGORITHM : baileyMul_fast<16,1,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,e;

   f16 = 0.0;
   f15 = ddg_two_prod(a_lolololo,b,&e);
   f16 += e;
   f14 = ddg_two_prod(a_hilololo,b,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f13 = ddg_two_prod(a_lohilolo,b,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f12 = ddg_two_prod(a_hihilolo,b,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f11 = ddg_two_prod(a_lolohilo,b,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f10 = ddg_two_prod(a_hilohilo,b,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f9 = ddg_two_prod(a_lohihilo,b,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f8 = ddg_two_prod(a_hihihilo,b,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f7 = ddg_two_prod(a_lololohi,b,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f6 = ddg_two_prod(a_hilolohi,b,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f5 = ddg_two_prod(a_lohilohi,b,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f4 = ddg_two_prod(a_hihilohi,b,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f3 = ddg_two_prod(a_lolohihi,b,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f2 = ddg_two_prod(a_hilohihi,b,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f1 = ddg_two_prod(a_lohihihi,b,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;
   f0 = ddg_two_prod(a_hihihihi,b,&e);
   f1 = ddg_two_sum(f1,e,&e);
   f2 = ddg_two_sum(f2,e,&e);
   f3 = ddg_two_sum(f3,e,&e);
   f4 = ddg_two_sum(f4,e,&e);
   f5 = ddg_two_sum(f5,e,&e);
   f6 = ddg_two_sum(f6,e,&e);
   f7 = ddg_two_sum(f7,e,&e);
   f8 = ddg_two_sum(f8,e,&e);
   f9 = ddg_two_sum(f9,e,&e);
   f10 = ddg_two_sum(f10,e,&e);
   f11 = ddg_two_sum(f11,e,&e);
   f12 = ddg_two_sum(f12,e,&e);
   f13 = ddg_two_sum(f13,e,&e);
   f14 = ddg_two_sum(f14,e,&e);
   f15 = ddg_two_sum(f15,e,&e);
   f16 += e;

   hdg_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                   c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
                   c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
                   c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
                   c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
}

__device__ __forceinline__ void hdg_mlt_d
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b )
{
   double c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi;
   double c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi;
   double c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo;
   double c_hihilolo,c_lohilolo,c_hilololo,c_lolololo;

   hdg_mul_hd_d(*a_hihihihi,*a_lohihihi,*a_hilohihi,*a_lolohihi,
                *a_hihilohi,*a_lohilohi,*a_hilolohi,*a_lololohi,
                *a_hihihilo,*a_lohihilo,*a_hilohilo,*a_lolohilo,
                *a_hihilolo,*a_lohilolo,*a_hilololo,*a_lolololo,b,
                &c_hihihihi,&c_lohihihi,&c_hilohihi,&c_lolohihi,
                &c_hihilohi,&c_lohilohi,&c_hilolohi,&c_lololohi,
                &c_hihihilo,&c_lohihilo,&c_hilohilo,&c_lolohilo,
                &c_hihilolo,&c_lohilolo,&c_hilololo,&c_lolololo);

   *a_hihihihi = c_hihihihi;
   *a_lohihihi = c_lohihihi;
   *a_hilohihi = c_hilohihi;
   *a_lolohihi = c_lolohihi;
   *a_hihilohi = c_hihilohi;
   *a_lohilohi = c_lohilohi;
   *a_hilolohi = c_hilolohi;
   *a_lololohi = c_lololohi;
   *a_hihihilo = c_hihihilo;
   *a_lohihilo = c_lohihilo;
   *a_hilohilo = c_hilohilo;
   *a_lolohilo = c_lolohilo;
   *a_hihilolo = c_hihilolo;
   *a_lohilolo = c_lohilolo;
   *a_hilololo = c_hilololo;
   *a_lolololo = c_lolololo;
}

__device__ __forceinline__ void hdg_div 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
   double q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16;
   double acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi;
   double acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi;
   double acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo;
   double acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo;

   q0 = a_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q0,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
           a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
           a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
           a_hihilolo,a_lohilolo,a_hilololo,a_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q1 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q1,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q2 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q2,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q3 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q3,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q4 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q4,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q5 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q5,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q6 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q6,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q7 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q7,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q8 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q8,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q9 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q9,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q10 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q10,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q11 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q11,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q12 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q12,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q13 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q13,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q14 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q14,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q15 = *c_hihihihi/b_hihihihi;
   hdg_mul_hd_d(b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo,q15,
                &acc_hihihihi,&acc_lohihihi,&acc_hilohihi,&acc_lolohihi,
                &acc_hihilohi,&acc_lohilohi,&acc_hilolohi,&acc_lololohi,
                &acc_hihihilo,&acc_lohihilo,&acc_hilohilo,&acc_lolohilo,
                &acc_hihilolo,&acc_lohilolo,&acc_hilololo,&acc_lolololo);
   hdg_sub(*c_hihihihi,*c_lohihihi,*c_hilohihi,*c_lolohihi,
           *c_hihilohi,*c_lohilohi,*c_hilolohi,*c_lololohi,
           *c_hihihilo,*c_lohihilo,*c_hilohilo,*c_lolohilo,
           *c_hihilolo,*c_lohilolo,*c_hilololo,*c_lolololo,
           acc_hihihihi,acc_lohihihi,acc_hilohihi,acc_lolohihi,
           acc_hihilohi,acc_lohilohi,acc_hilolohi,acc_lololohi,
           acc_hihihilo,acc_lohihilo,acc_hilohilo,acc_lolohilo,
           acc_hihilolo,acc_lohilolo,acc_hilololo,acc_lolololo,
           c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
           c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
           c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
           c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);

   q16 = *c_hihihihi/b_hihihihi;

   hdg_fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,
                   c_hihihihi,c_lohihihi,c_hilohihi,c_lolohihi,
                   c_hihilohi,c_lohilohi,c_hilolohi,c_lololohi,
                   c_hihihilo,c_lohihilo,c_hilohilo,c_lolohilo,
                   c_hihilolo,c_lohilolo,c_hilololo,c_lolololo);
}

/************************* square root *********************************/

__device__ __forceinline__ void hdg_sqrt 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
   double z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi;
   double z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi;
   double z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo;
   double z_hihilolo,z_lohilolo,z_hilololo,z_lolololo;

   odg_sqrt(a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
            a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
            b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
            b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi);
   hdg_sqr(*b_hihihihi,*b_lohihihi,*b_hilohihi,*b_lolohihi,
           *b_hihilohi,*b_lohilohi,*b_hilolohi,*b_lololohi,
           0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
           &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
           &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
           &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
           &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
   hdg_inc(&z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
           &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
           &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
           &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo,
           a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi,
           a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi,
           a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo,
           a_hihilolo,a_lohilolo,a_hilololo,a_lolololo);
   hdg_div(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
           z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
           z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
           z_hihilolo,z_lohilolo,z_hilololo,z_lolololo,
           *b_hihihihi,*b_lohihihi,*b_hilohihi,*b_lolohihi,
           *b_hihilohi,*b_lohilohi,*b_hilolohi,*b_lololohi,
           0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
           &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
           &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
           &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
           &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
   hdg_mul_pwr2(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
                z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
                z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
                z_hihilolo,z_lohilolo,z_hilololo,z_lolololo,0.5,
                b_hihihihi,b_lohihihi,b_hilohihi,b_lolohihi,
                b_hihilohi,b_lohilohi,b_hilolohi,b_lololohi,
                b_hihihilo,b_lohihilo,b_hilohilo,b_lolohilo,
                b_hihilolo,b_lohilolo,b_hilololo,b_lolololo);
}
