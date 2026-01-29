/* Collection of functions to split doubles. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "splitting_doubles.h"

using namespace std;

void write_52bits ( int k, uint64 nbr )
{
   if(k > 0)
   {
      write_52bits(k-1,nbr/2);
      cout << nbr % 2;
      if(k % 4 == 0) cout << " ";
   }
}

void write_52double ( double nbr )
{
   int exponent;
   double fraction = frexp(nbr, &exponent );
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;

   write_52bits(52, int64fac);
   cout << " " << exponent << endl;
}

uint64 last_bits ( int k, uint64 nbr )
{
   uint64 mask = 1;
   uint64 pwr = 2;

   for(int i=1; i<k; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   return (nbr & mask); // bit masking
}

void quarter_bits
 ( uint64 nbr, uint64 *b0, uint64 *b1, uint64 *b2, uint64 *b3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.quarter_bits ..." << endl;

   uint64 mask = 1;
   uint64 pwr = 2;
   uint64 rest = nbr;

   for(int i=1; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
   {
      cout << " m : "; write_52bits(52, mask); cout << endl;
   }
   *b3 = (nbr & mask); // last 13 bits
   rest = rest - *b3;
   
   if(vrblvl > 0)
   {
      cout << hex << " h : " << nbr << endl; cout << dec;
      cout << " b : "; write_52bits(52, nbr); cout << endl;
      cout << "b3 : "; write_52bits(52, *b3); cout << endl;
   }
   for(int i=0; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
   {
      cout << " m : "; write_52bits(52, mask); cout << endl;
   }
   *b2 = (rest & mask); // next 13 bits
   rest = rest - *b2;

   if(vrblvl > 0)
   {
      cout << "b2 : "; write_52bits(52,*b2); cout << endl;
   }
   for(int i=0; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
   {
      cout << " m : "; write_52bits(52, mask); cout << endl;
   }
   *b1 = (rest & mask); // next 13 bits

   if(vrblvl > 0)
   {
      cout << "b1 : "; write_52bits(52,*b1); cout << endl;
   }
   *b0 = nbr - *b1 - *b2 - *b3;

   if(vrblvl > 0)
   {
      cout << "b0 : "; write_52bits(52, *b0); cout << endl;
      rest = nbr - *b0 - *b1 - *b2 - *b3;
      cout << " r : "; write_52bits(52, rest);
      if(rest == 0)
         cout << " OKAY" << endl;
      else
         cout << " BUG!" << endl;
   }
}

double first_half ( double x, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.first_half ..." << endl;

   int exponent;
   double fraction = frexp(x, &exponent );
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;
   uint64 second_part = last_bits(26, int64fac);
   uint64 first_part = int64fac - second_part;
   double result = ldexp(first_part, exponent-52);

   if(vrblvl > 0)
   {
      cout << "exponent : " << exponent << endl;
      cout << "fraction : " << int64fac << endl;
      cout << hex << " h : " << int64fac << endl; cout << dec;
      cout << " b : "; write_52bits(52, int64fac); cout << endl;
      cout << "b0 : "; write_52bits(52, first_part); cout << endl;
      cout << "b1 : "; write_52bits(52, second_part); cout << endl;
   }
   return result;
}

void half_split ( double x, double *x0, double *x1, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.half_split ..." << endl;

   *x0 = first_half(x, vrblvl);
   *x1 = x - *x0;
}

int leading_zeros ( uint64 nbr, int idxpwr, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.leading_zeros ..." << endl;

   int result = 0;
   int idx = idxpwr;
   const uint64 two = 2;
   uint64 threshold = two << idx; // 2**idx

   if(vrblvl > 0)
      cout << "threshold : " << threshold << ", nbr : " << nbr << endl;

   while(nbr < threshold)
   {
      result = result + 1;
      threshold = threshold/2;
      if(vrblvl > 0)
         cout << "threshold : " << threshold << ", cnt : " << result << endl;
   }
   
   return result;
}

void quarter_split
 ( double x, double *x0, double *x1, double *x2, double *x3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.quarter_split ..." << endl;

   int exponent;
   double fraction = frexp(x, &exponent);
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;
   uint64 f0,f1,f2,f3;
   double xf0,xf1,xf2,xf3;
   int cnt;

   if(vrblvl > 0) cout << "exponent : " << exponent << endl;

   quarter_bits(int64fac, &f0, &f1, &f2, &f3, vrblvl);

   xf0 = (double) f0;
   *x0 = ldexp(xf0, exponent - 52);

   cnt = (f1 == 0) ? 0 : leading_zeros(f1, 37, vrblvl);
   if(vrblvl > 0) cout << "#leading zeros in f1 : " << cnt << endl;
   xf1 = (double) f1;
   // *x1 = ldexp(xf1, exponent - 52 - cnt); // does not work ...
   if(cnt == 0)
      *x1 = ldexp(xf1, exponent - 52);
   else
   {
      for(int i=0; i<cnt; i++) f1 = 2*f1;
      xf1 = (double) f1;
      *x1 = ldexp(xf1, exponent - 52);
      for(int i=0; i<cnt; i++) *x1 = *x1/2.0;
   }
   cnt = (f2 == 0) ? 0 : leading_zeros(f2, 24, vrblvl);
   if(vrblvl > 0) cout << "#leading zeros in f2 : " << cnt << endl;
   xf2 = (double) f2;
   if(cnt == 0)
      *x2 = ldexp(xf2, exponent - 52);
   else
   {
      for(int i=0; i<cnt; i++) f2 = 2*f2;
      xf2 = (double) f2;
      *x2 = ldexp(xf2, exponent - 52);
      for(int i=0; i<cnt; i++) *x2 = *x2/2.0;
   }
   if(vrblvl > 0) // just for testing purposes
   {
      cnt = (f3 == 0) ? 0 : leading_zeros(f3, 11, vrblvl);
      cout << "#leading zeros in f3 : " << cnt << endl;
   }
   *x3 = x - *x0 - *x1 - *x2;
}

bool is_quarter_balanced ( double x, double y, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.is_quarter_balanced ..." << endl;

   int xe,ye;
   double xf = frexp(x, &xe);
   double yf = frexp(y, &ye);
   int dxy = xe - ye;

   if(vrblvl > 0)
   {
      cout << "x : " << x << " has exponent " << xe << endl;
      cout << "y : " << y << " has exponent " << ye << endl;
      cout << "difference of exponents : " << dxy;
   }
   bool result = not (dxy > 13);
   if(vrblvl > 0)
   {
      if(result)
         cout << " balanced" << endl;
      else
         cout << " unbalanced" << endl;
   }
   return result;
}

void quarter_balance ( double *x, double *y, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.quarter_balance ..." << endl;

   int exn;
   double xf = frexp(*x, &exn);
   double bit = ldexp(1.0, exn - 13);

   if(vrblvl > 0)
   {
      cout << "b x : "; write_52double(*x);
      cout << "bit : "; write_52double(bit);
      cout << "b y : "; write_52double(*y);
   }
   *x = *x - bit;
   *y = *y + bit;

   if(vrblvl > 0)
   {
      cout << "b x : "; write_52double(*x);
      cout << "b y : "; write_52double(*y);
   }
}

void balance_quarters
 ( double *x0, double *x1, double *x2, double *x3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in splitting_doubles.balance_quarters ..." << endl;

   if(not is_quarter_balanced(*x0, *x1, vrblvl-1))
      quarter_balance(x0, x1, vrblvl-1);
   if(not is_quarter_balanced(*x1, *x2, vrblvl-1))
      quarter_balance(x1, x2, vrblvl-1);
   if(not is_quarter_balanced(*x2, *x3, vrblvl-1))
      quarter_balance(x2, x3, vrblvl-1);
}
