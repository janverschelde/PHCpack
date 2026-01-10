/* tests the removal of the last bits from a 64-bit double */

#include <stdio.h>
#include <math.h>

void write_52bits ( int* bits );
/*
 * Given in bits an array of 52 bits,
 * writes the bits in groups of four. */

unsigned long long int value_52bits ( int* bits );
/*
 * Given in bits an array of 52 bits,
 * returns its value as an unsigned 64-bit integer. */

void expand_52bits ( int* bits, unsigned long long int nbr );
/*
 * Given in bits space for 52 bits,
 * fills bits with the binary expansion of nbr. */

int main ( int argc, char *argv[] )
{
   const double x = M_PI;
   const int last = 26;  // remove the last bits from x
   int thebits[52];      // stores the bits of x

   printf("x = %.14le\n", x);

   int exponent;
   double fraction = frexp(x, &exponent );
   double shifted = ldexp(fraction, 52);
   unsigned long long int int64fac = (unsigned long long int) shifted;

   printf("fraction as a double : %.14e\n", fraction);
   printf("the exponent : %d\n", exponent);
   printf("fraction as 64-bit int : %lld\n", int64fac);
   printf("fraction in hex format : %llx\n", int64fac);

   double y = ldexp(fraction, exponent);

   printf("y = %.14le\n", y);

   expand_52bits(thebits, int64fac);
   printf("all bits :"); write_52bits(thebits); printf("\n");
   for(size_t i=0; i<last; i++) thebits[51-i] = 0;
   printf("all bits :"); write_52bits(thebits); printf("\n");

   unsigned long long int truncated = value_52bits(thebits);
   shifted = ldexp(truncated, -52);
   printf("truncated fraction : %.14e\n", shifted);
   int64fac = (unsigned long long int) ldexp(shifted, 52);
   printf("truncated fraction as 64-bit int : %lld\n", int64fac);
   printf("truncated fraction in hex format : %llx\n", int64fac);

   double z = ldexp(shifted, exponent);

   printf("z = %.14le\n", z); 
   printf("x = %.14le\n", x); 
   printf("error : %.2e\n", fabs(x - z));

   return 0;
}

void write_52bits ( int* bits )
{
   for(size_t i=0; i<52; i++)
   {
      if(i % 4 == 0) printf(" ");
      printf("%d", bits[i]);
   }
}

unsigned long long int value_52bits ( int* bits )
{
   unsigned long long int value = 0;

   for(size_t i=0; i<52; i++) value = 2*value + bits[i];

   return value;
}

void expand_52bits ( int* bits, unsigned long long int nbr )
{
   unsigned long long int temp = nbr;

   for(size_t i=0; i<52; i++)
   {
      bits[51-i] = temp % 2;
      temp = temp/2;
   }
}
