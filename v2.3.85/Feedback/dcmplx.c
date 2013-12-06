/* file "dcmplx.c" contains definitions of prototypes in "dcmplx.h" */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcmplx.h"



dcmplx create1 ( double r )
{
   dcmplx z;

   z.re = r;
   z.im = 0.0l;

   return z;
}

dcmplx create2 ( double r, double i )
{
   dcmplx z;

   z.re = r;
   z.im = i;

   return z;
}

dcmplx polar ( double r, double a )
{
   dcmplx z;

   z.re = r*cos(a);
   z.im = r*sin(a);

   return z;
}

dcmplx random_dcmplx1 ( void )
{
   dcmplx z;
   double angle = rand();

   z = polar(1,angle);

   return z;
}

void read_dcmplx ( dcmplx *z)
{
   double r,i;

   scanf("%lf %lf", &r, &i);

   z->re = r;
   z->im = i;
}

void read_dcmplx0 (dcmplx *z, FILE *ifp)
{
   double r; 
  
   fscanf(ifp, "%lf", &r);
   z->re = r;
   z->im = 0.0;
}

void read_dcmplx1 (dcmplx *z, FILE *ifp )
{
   double r,i;

   fscanf(ifp, "%lf %lf", &r, &i);

   z->re = r;
   z->im = i;
}
void write_dcmplx ( dcmplx z )
{
   printf("%.15le  %.15le*i  ", z.re, z.im);
}

void write_dcmplx1 ( dcmplx z, FILE *ofp )
{
   if(z.re >= 0) fprintf(ofp, " ");
   fprintf(ofp, "%.15le", z.re);
   if(z.im >= 0) fprintf(ofp, "+");
   fprintf(ofp, "%.15le*i", z.im);
}

void writeln_dcmplx ( dcmplx z )
{
   printf("%.15le  %.15le*i\n", z.re, z.im);
}

void writeln_dcmplx1 ( dcmplx z, FILE *ofp )
{
   if(z.re >= 0) fprintf(ofp, " ");
   fprintf(ofp, "%.15le", z.re);
   if(z.im >= 0) fprintf(ofp, "+");
   fprintf(ofp, "%.15le*i\n", z.im);
}

double dabs ( double x )
{
   if (x >= 0.0)
      return x;
   else
      return -x;
}

double dcabs ( dcmplx z )
{
   return (dabs(z.re) + dabs(z.im));
}

int equal_dcmplx ( dcmplx z1, dcmplx z2, double tol )
{
   if (dabs(z1.re - z2.re) > tol)
      return 0;
   else if (dabs(z1.im - z2.im) > tol)
      return 0;
   else
      return 1;
}

double modulus ( dcmplx z ) 
{
   double absre = dabs(z.re);
   double absim = dabs(z.im);

   if (absre >= absim)
      if (z.im == 0.0)
         return absre;
      else 
         return absre*sqrt(1.0+pow(z.im/z.re,2.0));
   else
      if (z.re == 0.0)
         return absim;
      else 
         return absim*sqrt(1.0+pow(z.re/z.im,2.0));
      
}

dcmplx conjugate ( dcmplx z )
{
   dcmplx cz;

   cz.re = z.re;
   cz.im = -z.im;

   return cz;
}

dcmplx min_dcmplx ( dcmplx z1 )
{
   dcmplx z;

   z.re = -z1.re;
   z.im = -z1.im;

   return z;
}

dcmplx add_dcmplx ( dcmplx z1, dcmplx z2 )
{
   dcmplx z;

   z.re = z1.re + z2.re;
   z.im = z1.im + z2.im;

   return z;
}

dcmplx sub_dcmplx ( dcmplx z1, dcmplx z2 )
{
   dcmplx z;

   z.re = z1.re - z2.re;
   z.im = z1.im - z2.im;

   return z;
}

dcmplx mul_dcmplx ( dcmplx z1, dcmplx z2 )
{
   dcmplx z;

   z.re = z1.re*z2.re - z1.im*z2.im;
   z.im = z1.im*z2.re + z1.re*z2.im;

   return z;
}

dcmplx div_dcmplx ( dcmplx z1, dcmplx z2 )
{
   dcmplx z;
   double absz2re,absz2im,tmp;

   if (z2.im == 0.0)
   {
      z.re = z1.re/z2.re;
      z.im = z1.im/z2.re;
   }
   else if (z2.re == 0.0)
   {   
      z.re = z1.im/z2.im;
      z.im = -z1.re/z2.im;
   }
   else
   {
      absz2re = dabs(z2.re);
      absz2im = dabs(z2.im);
      if (absz2re < absz2im)
      {
         tmp = z2.re/z2.im;
         z.re = z1.re*tmp; z.re += z1.im;
         z.im = z1.im*tmp; z.im -= z1.re;
         tmp *= z2.re; tmp += z2.im;
         z.re /= tmp; 
         z.im /= tmp; 
      }
      else if (absz2re > absz2im)
      {
         tmp = z2.im/z2.re;
         z.re = z1.im*tmp; z.re += z1.re;
         z.im = z1.re*tmp; z.im -= z1.im; z.im = -z.im;
         tmp *= z2.im; tmp += z2.re;
         z.re /= tmp; 
         z.im /= tmp;
      }
      else if (z2.re == z2.im)
      {
         tmp = 2*z2.re;
         z.re = z1.re + z1.im; z.re /= tmp;
         z.im = z1.im - z1.re; z.im /= tmp;
      }
      else
      {
         tmp = 2*z2.re;
         z.re = z1.re - z1.im; z.re /= tmp;
         z.im = z1.im + z1.re; z.im /= tmp;
      }
   }

   return z;
}

dcmplx add_double ( dcmplx z1, double z2 )
{
   dcmplx sum;

   sum.re = z1.re + z2;
   sum.im = z1.im;

   return sum;
}

dcmplx sub_double ( dcmplx z1, double z2 )
{
   dcmplx dif;

   dif.re = z1.re - z2;
   dif.im = z1.im;

   return dif;
}

dcmplx mul_double ( dcmplx z1, double z2 )
{
   dcmplx prod;

   prod.re = z1.re*z2;
   prod.im = z1.im*z2;

   return prod;
}

dcmplx div_double ( dcmplx z1, double z2 )
{
   dcmplx quot;

   quot.re = z1.re/z2;
   quot.im = z1.im/z2;

   return quot;
}

void swap(dcmplx* a, dcmplx* b)
{
  dcmplx tmp;
  tmp=*a;
  *a=*b;
  *b=tmp;
}
