#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dcmplx.h" 
#include "poly_gcd.h"
#include "poly_dcmplx.h"
#include "poly_smith.h"

void manual_test ( int n, int m );

int main(void)
{
  int n, m;

  printf("Please give the degree of the first polynomial:\n");
  scanf("%d", &n);
  
  printf("Please give the degree of the second polynomial:\n");
  scanf("%d", &m);

  manual_test(n, m);
  return 0;
}

void manual_test ( int n, int m )
{
  POLY a, b, c;
  int i;

  a.d=n;
  b.d=m;
  a.p=(dcmplx*) calloc(n+1, sizeof(dcmplx));
  b.p=(dcmplx*) calloc(m+1, sizeof(dcmplx));
  printf("please input the polynomial a\n"); 
  for(i=0; i<=n; i++)
     read_dcmplx( &(a.p[i]));
  printf("Polynomial a is:\n");
  Print_Poly(a.d, a.p);
  printf("please input the polynomial b\n");
  for(i=0; i<=m; i++)
     read_dcmplx(&(b.p[i]));
  printf("Polynomial b is:\n");
  Print_Poly(b.d, b.p);
  c=mul_poly1(a, b);
  printf("a*b = ");
  Print_Poly(c.d, c.p);

  free(a.p);
  free(b.p);
}
