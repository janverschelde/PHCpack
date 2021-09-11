// The file octo_double_functions.cpp defines the code for the functions
// specified in octo_double_functions.h

#include <cmath>
#include <iostream>
#include <iomanip>
#include "double_double_functions.h"
#include "quad_double_functions.h"
#include "octo_double_functions.h"

/************************** renormalizations **************************/

void odf_renorm8
 ( double f0, double f1, double f2, double f3, double f4, double f5,
   double f6, double f7, double f8, double *pr, double *r0, double *r1,
   double *r2, double *r3, double *r4, double *r5, double *r6, double *r7 )
{
   int ptr;

   if(f1 == 0.0)
   {
      *pr = f0;
      ptr = 0;
      *r0 = ddf_quick_two_sum(*pr,f2,pr);
   }
   else
   {
      *r0 = f0;
      *pr = f1;
      ptr = 1;
      *r1 = ddf_quick_two_sum(*pr,f2,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f3,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f3,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f4,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f4,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f5,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddf_quick_two_sum(*pr,f5,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddf_quick_two_sum(*pr,f6,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddf_quick_two_sum(*pr,f6,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddf_quick_two_sum(*pr,f7,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddf_quick_two_sum(*pr,f7,pr);
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
      *r0 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 1)
   {
      *r1 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 2)
   {
      *r2 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 3)
   {
      *r3 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 4)
   {
      *r4 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 5)
   {
      *r5 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 6)
   {
      *r6 = ddf_quick_two_sum(*pr,f8,pr);
   }
   else if(ptr == 7)
   {
      *r7 = ddf_quick_two_sum(*pr,f8,pr);
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
   if((ptr < 8) && (*pr != 0.0))
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
      ptr = ptr + 1;
   }
   if(ptr < 1)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0; *r3 = 0.0; *r2 = 0.0;
      *r1 = 0.0; *r0 = 0.0;
   }
   else if(ptr < 2)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0; *r3 = 0.0; *r2 = 0.0;
      *r1 = 0.0;
   }
   else if(ptr < 3)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0; *r3 = 0.0; *r2 = 0.0;
   }
   else if(ptr < 4)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0; *r3 = 0.0;
   }
   else if(ptr < 5)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0; *r4 = 0.0;
   }
   else if(ptr < 6)
   {
      *r7 = 0.0; *r6 = 0.0; *r5 = 0.0;
   }
   else if(ptr < 7)
   {
      *r7 = 0.0; *r6 = 0.0;
   }
   else if(ptr < 8)
   {
      *r7 = 0.0;
   }
}

void odf_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double *r0, double *r1, double *r2,
   double *r3, double *r4, double *r5, double *r6, double *r7 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,pr;

   pr = ddf_quick_two_sum(x7,x8,&f8);
   pr = ddf_quick_two_sum(x6,pr,&f7);
   pr = ddf_quick_two_sum(x5,pr,&f6);
   pr = ddf_quick_two_sum(x4,pr,&f5);
   pr = ddf_quick_two_sum(x3,pr,&f4);
   pr = ddf_quick_two_sum(x2,pr,&f3);
   pr = ddf_quick_two_sum(x1,pr,&f2);
   f0 = ddf_quick_two_sum(x0,pr,&f1);

   odf_renorm8(f0,f1,f2,f3,f4,f5,f6,f7,f8,&pr,r0,r1,r2,r3,r4,r5,r6,r7);
}

void odf_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double y, double *r0, double *r1, double *r2,
   double *r3, double *r4, double *r5, double *r6, double *r7 )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,pr;

   pr = ddf_two_sum(x7,y,&f8);
   pr = ddf_two_sum(x6,pr,&f7);
   pr = ddf_two_sum(x5,pr,&f6);
   pr = ddf_two_sum(x4,pr,&f5);
   pr = ddf_two_sum(x3,pr,&f4);
   pr = ddf_two_sum(x2,pr,&f3);
   pr = ddf_two_sum(x1,pr,&f2);
   f0 = ddf_two_sum(x0,pr,&f1);

   odf_renorm8(f0,f1,f2,f3,f4,f5,f6,f7,f8,&pr,r0,r1,r2,r3,r4,r5,r6,r7);
}

/************************ copy and abs *******************************/

void odf_copy
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double *b_hihihi, double *b_lohihi, double *b_hilohi, double *b_lolohi,
   double *b_hihilo, double *b_lohilo, double *b_hilolo, double *b_lololo )
{
   *b_hihihi = a_hihihi;
   *b_lohihi = a_lohihi;
   *b_hilohi = a_hilohi;
   *b_lolohi = a_lolohi;
   *b_hihilo = a_hihilo;
   *b_lohilo = a_lohilo;
   *b_hilolo = a_hilolo;
   *b_lololo = a_lololo;
}

void odf_abs
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double *b_hihihi, double *b_lohihi, double *b_hilohi, double *b_lolohi,
   double *b_hihilo, double *b_lohilo, double *b_hilolo, double *b_lololo )
{
   if(a_hihihi < 0.0)
   {
      *b_hihihi = -a_hihihi;
      *b_lohihi = -a_lohihi;
      *b_hilohi = -a_hilohi;
      *b_lolohi = -a_lolohi;
      *b_hihilo = -a_hihilo;
      *b_lohilo = -a_lohilo;
      *b_hilolo = -a_hilolo;
      *b_lololo = -a_lololo;
   }
   else
   {
      *b_hihihi = a_hihihi;
      *b_lohihi = a_lohihi;
      *b_hilohi = a_hilohi;
      *b_lolohi = a_lolohi;
      *b_hihilo = a_hihilo;
      *b_lohilo = a_lohilo;
      *b_hilolo = a_hilolo;
      *b_lololo = a_lololo;
   }
}

/****************** additions and substractions ************************/

void odf_add
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   // ALGORITHM : baileyAddf_fast<8,8,8> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,e;

   f8 = 0.0;
   f7 = ddf_two_sum(a_lololo,b_lololo,&e);
   f8 += e;
   f6 = ddf_two_sum(a_hilolo,b_hilolo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(a_lohilo,b_lohilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(a_hihilo,b_hihilo,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(a_lolohi,b_lolohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(a_hilohi,b_hilohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_sum(a_lohihi,b_lohihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f0 = ddf_two_sum(a_hihihi,b_hihihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;

   odf_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                   c_hihihi,c_lohihi,c_hilohi,c_lolohi,
                   c_hihilo,c_lohilo,c_hilolo,c_lololo);
}

void odf_inc
 ( double *a_hihihi, double *a_lohihi, double *a_hilohi, double *a_lolohi,
   double *a_hihilo, double *a_lohilo, double *a_hilolo, double *a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,e;

   f8 = 0.0;
   f7 = ddf_two_sum(*a_lololo,b_lololo,&e);
   f8 += e;
   f6 = ddf_two_sum(*a_hilolo,b_hilolo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(*a_lohilo,b_lohilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(*a_hihilo,b_hihilo,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(*a_lolohi,b_lolohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(*a_hilohi,b_hilohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_sum(*a_lohihi,b_lohihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f0 = ddf_two_sum(*a_hihihi,b_hihihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;

   odf_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                   a_hihihi,a_lohihi,a_hilohi,a_lolohi,
                   a_hihilo,a_lohilo,a_hilolo,a_lololo);
}

void odf_inc_d
 ( double *a_hihihi, double *a_lohihi, double *a_hilohi, double *a_lolohi,
   double *a_hihilo, double *a_lohilo, double *a_hilolo, double *a_lololo,
   double b )
{
   odf_renorm_add1(*a_hihihi,*a_lohihi,*a_hilohi,*a_lolohi,
                   *a_hihilo,*a_lohilo,*a_hilolo,*a_lololo,b,
                   a_hihihi,a_lohihi,a_hilohi,a_lolohi,
                   a_hihilo,a_lohilo,a_hilolo,a_lololo);
}

void odf_dec
 ( double *a_hihihi, double *a_lohihi, double *a_hilohi, double *a_lolohi,
   double *a_hihilo, double *a_lohilo, double *a_hilolo, double *a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo )
{
   double mbhihihi = -b_hihihi;
   double mblohihi = -b_lohihi;
   double mbhilohi = -b_hilohi;
   double mblolohi = -b_lolohi;
   double mbhihilo = -b_hihilo;
   double mblohilo = -b_lohilo;
   double mbhilolo = -b_hilolo;
   double mblololo = -b_lololo;

   odf_inc(a_hihihi,a_lohihi,a_hilohi,a_lolohi,
           a_hihilo,a_lohilo,a_hilolo,a_lololo,
           mbhihihi,mblohihi,mbhilohi,mblolohi,
           mbhihilo,mblohilo,mbhilolo,mblololo);
}

void odf_minus
 ( double *a_hihihi, double *a_lohihi, double *a_hilohi, double *a_lolohi,
   double *a_hihilo, double *a_lohilo, double *a_hilolo, double *a_lololo )
{

   *a_hihihi = -(*a_hihihi);
   *a_lohihi = -(*a_lohihi);
   *a_hilohi = -(*a_hilohi);
   *a_lolohi = -(*a_lolohi);
   *a_hihilo = -(*a_hihilo);
   *a_lohilo = -(*a_lohilo);
   *a_hilolo = -(*a_hilolo);
   *a_lololo = -(*a_lololo);
}

void odf_sub
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   odf_copy(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
            b_hihilo,b_lohilo,b_hilolo,b_lololo,
            c_hihihi,c_lohihi,c_hilohi,c_lolohi,
            c_hihilo,c_lohilo,c_hilolo,c_lololo);
   odf_minus(c_hihihi,c_lohihi,c_hilohi,c_lolohi,
             c_hihilo,c_lohilo,c_hilolo,c_lololo);
   odf_inc(c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo,
           a_hihihi,a_lohihi,a_hilohi,a_lolohi,
           a_hihilo,a_lohilo,a_hilolo,a_lololo);
}

/***************** multiplications and division ********************/

void odf_mul_pwr2
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   *c_hihihi = a_hihihi*b;
   *c_lohihi = a_lohihi*b;
   *c_hilohi = a_hilohi*b;
   *c_lolohi = a_lolohi*b;
   *c_hihilo = a_hihilo*b;
   *c_lohilo = a_lohilo*b;
   *c_hilolo = a_hilolo*b;
   *c_lololo = a_lololo*b;
}

void odf_mul
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   // ALGORITHM : baileyMul_fast<8,8,8> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,p,e;

   f8 =  a_lohihi*b_lololo;
   f8 += a_hilohi*b_hilolo;
   f8 += a_lolohi*b_lohilo;
   f8 += a_hihilo*b_hihilo;
   f8 += a_lohilo*b_lolohi;
   f8 += a_hilolo*b_hilohi;
   f8 += a_lololo*b_lohihi;
   f7 = ddf_two_prod(a_hihihi,b_lololo,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_hilolo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_lohilo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,b_hihilo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,b_lolohi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,b_hilohi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hilolo,b_lohihi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lololo,b_hihihi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   f6 = ddf_two_prod(a_hihihi,b_hilolo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_lohilo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_hihilo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,b_lolohi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,b_hilohi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,b_lohihi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilolo,b_hihihi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_prod(a_hihihi,b_lohilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_hihilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_lolohi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,b_hilohi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,b_lohihi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,b_hihihi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_prod(a_hihihi,b_hihilo,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_lolohi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_hilohi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,b_lohihi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,b_hihihi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_prod(a_hihihi,b_lolohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_hilohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_lohihi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,b_hihihi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_prod(a_hihihi,b_hilohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_lohihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,b_hihihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_prod(a_hihihi,b_lohihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,b_hihihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f0 = ddf_two_prod(a_hihihi,b_hihihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;

   odf_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                   c_hihihi,c_lohihi,c_hilohi,c_lolohi,
                   c_hihilo,c_lohilo,c_hilolo,c_lololo);
}

void odf_sqr
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,p,e;

   f8 =  a_lohihi*a_lololo;
   f8 += a_hilohi*a_hilolo;
   f8 += a_lolohi*a_lohilo;
   f8 += a_hihilo*a_hihilo;
   f8 += a_lohilo*a_lolohi;
   f8 += a_hilolo*a_hilohi;
   f8 += a_lololo*a_lohihi;
   f7 = ddf_two_prod(a_hihihi,a_lololo,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_hilolo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_lohilo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,a_hihilo,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,a_lolohi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,a_hilohi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_hilolo,a_lohihi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   p = ddf_two_prod(a_lololo,a_hihihi,&e);
   f8 += e;
   f7 = ddf_two_sum(f7,p,&e);
   f8 += e;
   f6 = ddf_two_prod(a_hihihi,a_hilolo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_lohilo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_hihilo,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,a_lolohi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,a_hilohi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,a_lohihi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilolo,a_hihihi,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f6 = ddf_two_sum(f6,p,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_prod(a_hihihi,a_lohilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_hihilo,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_lolohi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,a_hilohi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,a_lohihi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohilo,a_hihihi,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_sum(f5,p,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_prod(a_hihihi,a_hihilo,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_lolohi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_hilohi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,a_lohihi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hihilo,a_hihihi,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_sum(f4,p,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_prod(a_hihihi,a_lolohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_hilohi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_lohihi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lolohi,a_hihihi,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_sum(f3,p,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_prod(a_hihihi,a_hilohi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_lohihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_hilohi,a_hihihi,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_sum(f2,p,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_prod(a_hihihi,a_lohihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   p = ddf_two_prod(a_lohihi,a_hihihi,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_sum(f1,p,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f0 = ddf_two_prod(a_hihihi,a_hihihi,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;

   odf_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                   c_hihihi,c_lohihi,c_hilohi,c_lolohi,
                   c_hihilo,c_lohilo,c_hilolo,c_lololo);
}

void odf_mul_od_d
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   // ALGORITHM : baileyMul_fast<8,1,8> generated by CAMPARY.

   double f0,f1,f2,f3,f4,f5,f6,f7,f8,e;

   f8 = 0.0;
   f7 = ddf_two_prod(a_lololo,b,&e);
   f8 += e;
   f6 = ddf_two_prod(a_hilolo,b,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f5 = ddf_two_prod(a_lohilo,b,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f4 = ddf_two_prod(a_hihilo,b,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f3 = ddf_two_prod(a_lolohi,b,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f2 = ddf_two_prod(a_hilohi,b,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f1 = ddf_two_prod(a_lohihi,b,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;
   f0 = ddf_two_prod(a_hihihi,b,&e);
   f1 = ddf_two_sum(f1,e,&e);
   f2 = ddf_two_sum(f2,e,&e);
   f3 = ddf_two_sum(f3,e,&e);
   f4 = ddf_two_sum(f4,e,&e);
   f5 = ddf_two_sum(f5,e,&e);
   f6 = ddf_two_sum(f6,e,&e);
   f7 = ddf_two_sum(f7,e,&e);
   f8 += e;

   odf_fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                   c_hihihi,c_lohihi,c_hilohi,c_lolohi,
                   c_hihilo,c_lohilo,c_hilolo,c_lololo);
}

void odf_mlt_d
 ( double *a_hihihi, double *a_lohihi, double *a_hilohi, double *a_lolohi,
   double *a_hihilo, double *a_lohilo, double *a_hilolo, double *a_lololo,
   double b )
{
   double c_hihihi,c_lohihi,c_hilohi,c_lolohi;
   double c_hihilo,c_lohilo,c_hilolo,c_lololo;

   odf_mul_od_d(*a_hihihi,*a_lohihi,*a_hilohi,*a_lolohi,
                *a_hihilo,*a_lohilo,*a_hilolo,*a_lololo,b,
                &c_hihihi,&c_lohihi,&c_hilohi,&c_lolohi,
                &c_hihilo,&c_lohilo,&c_hilolo,&c_lololo);

   *a_hihihi = c_hihihi;
   *a_lohihi = c_lohihi;
   *a_hilohi = c_hilohi;
   *a_lolohi = c_lolohi;
   *a_hihilo = c_hihilo;
   *a_lohilo = c_lohilo;
   *a_hilolo = c_hilolo;
   *a_lololo = c_lololo;
}

void odf_div
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double b_hihihi, double b_lohihi, double b_hilohi, double b_lolohi,
   double b_hihilo, double b_lohilo, double b_hilolo, double b_lololo,
   double *c_hihihi, double *c_lohihi, double *c_hilohi, double *c_lolohi,
   double *c_hihilo, double *c_lohilo, double *c_hilolo, double *c_lololo )
{
   double q0,q1,q2,q3,q4,q5,q6,q7,q8;
   double acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi;
   double acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo;

   q0 = a_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q0,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(a_hihihi,a_lohihi,a_hilohi,a_lolohi,
           a_hihilo,a_lohilo,a_hilolo,a_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q1 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q1,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q2 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q2,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q3 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q3,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q4 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q4,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q5 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q5,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q6 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q6,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q7 = *c_hihihi/b_hihihi;
   odf_mul_od_d(b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo,q7,
                &acc_hihihi,&acc_lohihi,&acc_hilohi,&acc_lolohi,
                &acc_hihilo,&acc_lohilo,&acc_hilolo,&acc_lololo);
   odf_sub(*c_hihihi,*c_lohihi,*c_hilohi,*c_lolohi,
           *c_hihilo,*c_lohilo,*c_hilolo,*c_lololo,
           acc_hihihi,acc_lohihi,acc_hilohi,acc_lolohi,
           acc_hihilo,acc_lohilo,acc_hilolo,acc_lololo,
           c_hihihi,c_lohihi,c_hilohi,c_lolohi,
           c_hihilo,c_lohilo,c_hilolo,c_lololo);

   q8 = *c_hihihi/b_hihihi;

   odf_fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,
                   c_hihihi,c_lohihi,c_hilohi,c_lolohi,
                   c_hihilo,c_lohilo,c_hilolo,c_lololo);
}

/***************************** square root *****************************/

void odf_sqrt
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo,
   double *b_hihihi, double *b_lohihi, double *b_hilohi, double *b_lolohi,
   double *b_hihilo, double *b_lohilo, double *b_hilolo, double *b_lololo )
{
   double z_hihihi,z_lohihi,z_hilohi,z_lolohi;
   double z_hihilo,z_lohilo,z_hilolo,z_lololo;

   qdf_sqrt(a_hihihi,a_lohihi,a_hilohi,a_lolohi,
            b_hihihi,b_lohihi,b_hilohi,b_lolohi);
   odf_sqr(*b_hihihi,*b_lohihi,*b_hilohi,*b_lolohi,0.0,0.0,0.0,0.0,
           &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
           &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
   odf_inc(&z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
           &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo,
           a_hihihi,a_lohihi,a_hilohi,a_lolohi,
           a_hihilo,a_lohilo,a_hilolo,a_lololo);
   odf_div(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
           z_hihilo,z_lohilo,z_hilolo,z_lololo,
           *b_hihihi,*b_lohihi,*b_hilohi,*b_lolohi,0.0,0.0,0.0,0.0,
           &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
           &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
   odf_mul_pwr2(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
                z_hihilo,z_lohilo,z_hilolo,z_lololo,0.5,
                b_hihihi,b_lohihi,b_hilohi,b_lolohi,
                b_hihilo,b_lohilo,b_hilolo,b_lololo);
}

/*************************** basic output ***************************/

void odf_write_doubles
 ( double a_hihihi, double a_lohihi, double a_hilohi, double a_lolohi,
   double a_hihilo, double a_lohilo, double a_hilolo, double a_lololo )
{
   std::cout << std::scientific << std::setprecision(16);
   std::cout << "  hihihi : " << a_hihihi;
   std::cout << "  lohihi : " << a_lohihi << std::endl;
   std::cout << "  hilohi : " << a_hilohi;
   std::cout << "  lolohi : " << a_lolohi << std::endl;
   std::cout << "  hihilo : " << a_hihilo;
   std::cout << "  lohilo : " << a_lohilo << std::endl;
   std::cout << "  hilolo : " << a_hilolo;
   std::cout << "  lololo : " << a_lololo << std::endl;
}
