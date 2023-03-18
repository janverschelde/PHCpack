/* file: hexa_double.c */

/* This file contains the corresponding C code for the functions
   with prototypes declared in the hexa_double.h file. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "double_double.h"
#include "hexa_double.h"

/************************* normalizations ************************/

void hd_renorm16
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
      *r0 = dd_quick_two_sum(*pr,f2,pr);
   }
   else
   {
      *r0 = f0;
      *pr = f1;
      ptr = 1;
      *r1 = dd_quick_two_sum(*pr,f2,pr);
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
      *r0 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f11,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f11,pr);
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
      *r0 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f12,pr);
   }
   else if(ptr == 11)
   {
      *r11 = dd_quick_two_sum(*pr,f12,pr);
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
      *r0 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 11)
   {
      *r11 = dd_quick_two_sum(*pr,f13,pr);
   }
   else if(ptr == 12)
   {
      *r12 = dd_quick_two_sum(*pr,f13,pr);
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
      *r0 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 11)
   {
      *r11 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 12)
   {
      *r12 = dd_quick_two_sum(*pr,f14,pr);
   }
   else if(ptr == 13)
   {
      *r13 = dd_quick_two_sum(*pr,f14,pr);
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
      *r0 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 11)
   {
      *r11 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 12)
   {
      *r12 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 13)
   {
      *r13 = dd_quick_two_sum(*pr,f15,pr);
   }
   else if(ptr == 14)
   {
      *r14 = dd_quick_two_sum(*pr,f15,pr);
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
      *r0 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 1)
   {
      *r1 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 2)
   {
      *r2 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 3)
   {
      *r3 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 4)
   {
      *r4 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 5)
   {
      *r5 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 6)
   {
      *r6 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 7)
   {
      *r7 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 8)
   {
      *r8 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 9)
   {
      *r9 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 10)
   {
      *r10 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 11)
   {
      *r11 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 12)
   {
      *r12 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 13)
   {
      *r13 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 14)
   {
      *r14 = dd_quick_two_sum(*pr,f16,pr);
   }
   else if(ptr == 15)
   {
      *r15 = dd_quick_two_sum(*pr,f16,pr);
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

void hd_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double x16,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr;

   pr = dd_quick_two_sum(x15,x16,&f16);
   pr = dd_quick_two_sum(x14,pr,&f15);
   pr = dd_quick_two_sum(x13,pr,&f14);
   pr = dd_quick_two_sum(x12,pr,&f13);
   pr = dd_quick_two_sum(x11,pr,&f12);
   pr = dd_quick_two_sum(x10,pr,&f11);
   pr = dd_quick_two_sum(x9,pr,&f10);
   pr = dd_quick_two_sum(x8,pr,&f9);
   pr = dd_quick_two_sum(x7,pr,&f8);
   pr = dd_quick_two_sum(x6,pr,&f7);
   pr = dd_quick_two_sum(x5,pr,&f6);
   pr = dd_quick_two_sum(x4,pr,&f5);
   pr = dd_quick_two_sum(x3,pr,&f4);
   pr = dd_quick_two_sum(x2,pr,&f3);
   pr = dd_quick_two_sum(x1,pr,&f2);
   f0 = dd_quick_two_sum(x0,pr,&f1);

   hd_renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,&pr,
               r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
}

void hd_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr;

   pr = dd_two_sum(x15,y,&f16);
   pr = dd_two_sum(x14,pr,&f15);
   pr = dd_two_sum(x13,pr,&f14);
   pr = dd_two_sum(x12,pr,&f13);
   pr = dd_two_sum(x11,pr,&f12);
   pr = dd_two_sum(x10,pr,&f11);
   pr = dd_two_sum(x9,pr,&f10);
   pr = dd_two_sum(x8,pr,&f9);
   pr = dd_two_sum(x7,pr,&f8);
   pr = dd_two_sum(x6,pr,&f7);
   pr = dd_two_sum(x5,pr,&f6);
   pr = dd_two_sum(x4,pr,&f5);
   pr = dd_two_sum(x3,pr,&f4);
   pr = dd_two_sum(x2,pr,&f3);
   pr = dd_two_sum(x1,pr,&f2);
   f0 = dd_two_sum(x0,pr,&f1);

   hd_renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,&pr,
               r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
}

/****************************** copy *****************************/

void hd_copy ( const double *a, double *b )
{
   b[0]  = a[0];  b[1]  = a[1];  b[2]  = a[2];  b[3]  = a[3];
   b[4]  = a[4];  b[5]  = a[5];  b[6]  = a[6];  b[7]  = a[7];
   b[8]  = a[8];  b[9]  = a[9];  b[10] = a[10]; b[11] = a[11];
   b[12] = a[12]; b[13] = a[13]; b[14] = a[14]; b[15] = a[15];
}

/******************* addition and subtraction *********************/

void hd_add ( const double *a, const double *b, double *c )
{
   // ALGORITHM : baileyAdd_fast<16,16,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,e;

   f16 = 0.0;
   f15 = dd_two_sum(a[15],b[15],&e);
   f16 += e;
   f14 = dd_two_sum(a[14],b[14],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(a[13],b[13],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(a[12],b[12],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(a[11],b[11],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(a[10],b[10],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(a[9],b[9],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(a[8],b[8],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(a[7],b[7],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(a[6],b[6],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(a[5],b[5],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_sum(a[4],b[4],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f3 = dd_two_sum(a[3],b[3],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f2 = dd_two_sum(a[2],b[2],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f1 = dd_two_sum(a[1],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
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
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   hd_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                  &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],
                  &c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15]);
}

void hd_add_hd_d ( const double *a, double b, double *c )
{
   hd_renorm_add1(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],
                  a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15],b,
                  &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],
                  &c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15]);
}

void hd_minus ( double *a )
{
   a[0]  = -a[0];  a[1]  = -a[1];  a[2]  = -a[2];  a[3]  = -a[3];
   a[4]  = -a[4];  a[5]  = -a[5];  a[6]  = -a[6];  a[7]  = -a[7];
   a[8]  = -a[8];  a[9]  = -a[9];  a[10] = -a[10]; a[11] = -a[11];
   a[12] = -a[12]; a[13] = -a[13]; a[14] = -a[14]; a[15] = -a[15];
}

void hd_sub ( const double *a, const double *b, double *c )
{
   double minb[16];

   hd_copy(b,minb);
   hd_minus(minb);
   hd_add(a,minb,c);
}

/**************  multiplications and division ***********************/

void hd_mul ( const double *a, const double *b, double *c )
{
   // ALGORITHM : baileyMul_fast<16,16,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,p,e;

   f16 =  a[1]*b[15];
   f16 += a[2]*b[14];
   f16 += a[3]*b[13];
   f16 += a[4]*b[12];
   f16 += a[5]*b[11];
   f16 += a[6]*b[10];
   f16 += a[7]*b[9];
   f16 += a[8]*b[8];
   f16 += a[9]*b[7];
   f16 += a[10]*b[6];
   f16 += a[11]*b[5];
   f16 += a[12]*b[4];
   f16 += a[13]*b[3];
   f16 += a[14]*b[2];
   f16 += a[15]*b[1];

   f15 = dd_two_prod(a[0],b[15],&e);
   f16 += e;
   p = dd_two_prod(a[1],b[14],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[13],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[12],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[11],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[10],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[9],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[8],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[7],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[6],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[5],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[11],b[4],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[12],b[3],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[13],b[2],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[14],b[1],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;
   p = dd_two_prod(a[15],b[0],&e);
   f16 += e;
   f15 = dd_two_sum(f15,p,&e);
   f16 += e;

   f14 = dd_two_prod(a[0],b[14],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[13],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[12],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[11],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[10],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[9],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[8],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[7],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[6],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[5],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[4],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[11],b[3],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[12],b[2],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[13],b[1],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[14],b[0],&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f14 = dd_two_sum(f14,p,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f13 = dd_two_prod(a[0],b[13],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[12],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[11],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[10],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[9],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[8],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[7],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[6],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[5],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[4],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[3],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[11],b[2],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[12],b[1],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[13],b[0],&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_sum(f13,p,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f12 = dd_two_prod(a[0],b[12],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[11],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[10],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[9],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[8],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[7],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[6],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[5],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[4],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[3],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[2],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[11],b[1],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[12],b[0],&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_sum(f12,p,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f11 = dd_two_prod(a[0],b[11],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[10],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[9],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[8],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[7],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[6],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[5],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[4],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[3],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[2],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[1],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[11],b[0],&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_sum(f11,p,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f10 = dd_two_prod(a[0],b[10],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[9],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[8],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[7],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[6],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[5],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[4],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[3],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[2],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[1],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[10],b[0],&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_sum(f10,p,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f9 = dd_two_prod(a[0],b[9],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[8],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[7],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[6],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[5],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[4],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[3],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[2],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[1],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[9],b[0],&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_sum(f9,p,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f8 = dd_two_prod(a[0],b[8],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[7],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[6],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[5],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[4],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[3],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[2],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[1],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[8],b[0],&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_sum(f8,p,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f7 = dd_two_prod(a[0],b[7],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[6],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[5],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[4],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[3],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[2],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[1],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[7],b[0],&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_sum(f7,p,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f6 = dd_two_prod(a[0],b[6],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[5],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[4],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[3],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[2],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[1],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[6],b[0],&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_sum(f6,p,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f5 = dd_two_prod(a[0],b[5],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[4],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[3],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[2],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[1],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[5],b[0],&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_sum(f5,p,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f4 = dd_two_prod(a[0],b[4],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[3],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[2],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[1],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[4],b[0],&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_sum(f4,p,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f3 = dd_two_prod(a[0],b[3],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[2],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[1],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[3],b[0],&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f3 = dd_two_sum(f3,p,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f2 = dd_two_prod(a[0],b[2],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[1],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[2],b[0],&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f2 = dd_two_sum(f2,p,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   f1 = dd_two_prod(a[0],b[1],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   p = dd_two_prod(a[1],b[0],&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f1 = dd_two_sum(f1,p,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

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
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   hd_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                  &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],
                  &c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15]);
}

void hd_mul_hd_d ( const double *a, double b, double *c )
{
   // ALGORITHM : baileyMul_fast<16,1,16> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,e;

   f16 = 0.0;
   f15 = dd_two_prod(a[15],b,&e);
   f16 += e;
   f14 = dd_two_prod(a[14],b,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f13 = dd_two_prod(a[13],b,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f12 = dd_two_prod(a[12],b,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f11 = dd_two_prod(a[11],b,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f10 = dd_two_prod(a[10],b,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f9 = dd_two_prod(a[9],b,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f8 = dd_two_prod(a[8],b,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f7 = dd_two_prod(a[7],b,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f6 = dd_two_prod(a[6],b,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f5 = dd_two_prod(a[5],b,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f4 = dd_two_prod(a[4],b,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f3 = dd_two_prod(a[3],b,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f2 = dd_two_prod(a[2],b,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
   f1 = dd_two_prod(a[1],b,&e);
   f2 = dd_two_sum(f2,e,&e);
   f3 = dd_two_sum(f3,e,&e);
   f4 = dd_two_sum(f4,e,&e);
   f5 = dd_two_sum(f5,e,&e);
   f6 = dd_two_sum(f6,e,&e);
   f7 = dd_two_sum(f7,e,&e);
   f8 = dd_two_sum(f8,e,&e);
   f9 = dd_two_sum(f9,e,&e);
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;
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
   f10 = dd_two_sum(f10,e,&e);
   f11 = dd_two_sum(f11,e,&e);
   f12 = dd_two_sum(f12,e,&e);
   f13 = dd_two_sum(f13,e,&e);
   f14 = dd_two_sum(f14,e,&e);
   f15 = dd_two_sum(f15,e,&e);
   f16 += e;

   hd_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                  &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],
                  &c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15]);
}

void hd_div ( const double *a, const double *b, double *c )
{
   double acc[16];
   double q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16;

   q0  = a[0]/b[0]; hd_mul_hd_d(b,q0,acc);  hd_sub(a,acc,c);
   q1  = c[0]/b[0]; hd_mul_hd_d(b,q1,acc);  hd_sub(c,acc,c);
   q2  = c[0]/b[0]; hd_mul_hd_d(b,q2,acc);  hd_sub(c,acc,c);
   q3  = c[0]/b[0]; hd_mul_hd_d(b,q3,acc);  hd_sub(c,acc,c);
   q4  = c[0]/b[0]; hd_mul_hd_d(b,q4,acc);  hd_sub(c,acc,c);
   q5  = c[0]/b[0]; hd_mul_hd_d(b,q5,acc);  hd_sub(c,acc,c);
   q6  = c[0]/b[0]; hd_mul_hd_d(b,q6,acc);  hd_sub(c,acc,c);
   q7  = c[0]/b[0]; hd_mul_hd_d(b,q7,acc);  hd_sub(c,acc,c);
   q8  = c[0]/b[0]; hd_mul_hd_d(b,q8,acc);  hd_sub(c,acc,c);
   q9  = c[0]/b[0]; hd_mul_hd_d(b,q9,acc);  hd_sub(c,acc,c);
   q10 = c[0]/b[0]; hd_mul_hd_d(b,q10,acc); hd_sub(c,acc,c);
   q11 = c[0]/b[0]; hd_mul_hd_d(b,q11,acc); hd_sub(c,acc,c);
   q12 = c[0]/b[0]; hd_mul_hd_d(b,q12,acc); hd_sub(c,acc,c);
   q13 = c[0]/b[0]; hd_mul_hd_d(b,q13,acc); hd_sub(c,acc,c);
   q14 = c[0]/b[0]; hd_mul_hd_d(b,q14,acc); hd_sub(c,acc,c);
   q15 = c[0]/b[0]; hd_mul_hd_d(b,q15,acc); hd_sub(c,acc,c);
   q16 = c[0]/b[0];

   hd_fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,
                  &c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],
                  &c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15]);
}

/******************* random number generator ***************************/

void hd_random ( double *x )
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
   const double eps10 = eps*eps9;
   const double eps11 = eps*eps10;
   const double eps12 = eps*eps11;
   const double eps13 = eps*eps12;
   const double eps14 = eps*eps13;
   const double eps15 = eps*eps14;
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
   const double r10 = ((double) rand())/RAND_MAX;
   const double r11 = ((double) rand())/RAND_MAX;
   const double r12 = ((double) rand())/RAND_MAX;
   const double r13 = ((double) rand())/RAND_MAX;
   const double r14 = ((double) rand())/RAND_MAX;
   const double r15 = ((double) rand())/RAND_MAX;

   x[0]  = r0;  x[1]  = 0.0; x[2]  = 0.0; x[3]  = 0.0;
   x[4]  = 0.0; x[5]  = 0.0; x[6]  = 0.0; x[7]  = 0.0;
   x[8]  = 0.0; x[9]  = 0.0; x[10] = 0.0; x[11] = 0.0;
   x[12] = 0.0; x[13] = 0.0; x[14] = 0.0; x[15] = 0.0;

   hd_add_hd_d(x,r1*eps,x);
   hd_add_hd_d(x,r2*eps2,x);
   hd_add_hd_d(x,r3*eps3,x);
   hd_add_hd_d(x,r4*eps4,x);
   hd_add_hd_d(x,r5*eps5,x);
   hd_add_hd_d(x,r6*eps6,x);
   hd_add_hd_d(x,r7*eps7,x);
   hd_add_hd_d(x,r8*eps8,x);
   hd_add_hd_d(x,r9*eps9,x);
   hd_add_hd_d(x,r10*eps10,x);
   hd_add_hd_d(x,r11*eps11,x);
   hd_add_hd_d(x,r12*eps12,x);
   hd_add_hd_d(x,r13*eps13,x);
   hd_add_hd_d(x,r14*eps14,x);
   hd_add_hd_d(x,r15*eps15,x);
}

/************************ basic output *********************************/

void hd_write_doubles ( const double *x )
{
   printf("  hihihihi = %21.14e",x[0]);
   printf("  lohihihi = %21.14e\n",x[1]);
   printf("  hilohihi = %21.14e",x[2]);
   printf("  lolohihi = %21.14e\n",x[3]);
   printf("  hihilohi = %21.14e",x[4]);
   printf("  lohilohi = %21.14e\n",x[5]);
   printf("  hilolohi = %21.14e",x[6]);
   printf("  lololohi = %21.14e\n",x[7]);
   printf("  hihihilo = %21.14e",x[8]);
   printf("  lohihilo = %21.14e\n",x[9]);
   printf("  hilohilo = %21.14e",x[10]);
   printf("  lolohilo = %21.14e\n",x[11]);
   printf("  hihilolo = %21.14e",x[12]);
   printf("  lohilolo = %21.14e\n",x[13]);
   printf("  hilololo = %21.14e",x[14]);
   printf("  lolololo = %21.14e\n",x[15]);
}
