#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dcmplx.h" 
#include "dc_roots.h"
#include "poly_gcd.h"
#include "poly_dcmplx.h"

#include "pieri_sols.h"

void random_test ( int n, int m, int l );
void manual_test ( int n, int m, int l );

int main(void)
{
  int k, n, m, l, i, t;

  printf("How many tests you want: \n");
  scanf("%d", &k);
  printf("Please give the degree of the first polynomial:\n");
  scanf("%d", &n);
  printf("Please give the degree of the second polynomial:\n");
  scanf("%d", &m);
  printf("How many common divisors you want in the test: \n");
  scanf("%d", &l);

  printf("Please choose what kind of test you want:\n");
  printf("1. random test\n");
  printf("2. interactive test\n");
  scanf("%d", &t);

  for(i=1; i<=k; i++)
  {
    printf("\n******* The No. %d test *********\n", i);
    if(t==1) random_test (n, m, l);  
    if(t==2) manual_test (n, m, l); 
  }

 return(0);
}



void manual_test ( int n, int m, int l )
{
  double eps = 1.0e-12;
  double tol = 1.0e-8;
  double mul_tol = 1.0e-3;
  dcmplx aa[n+1], bb[m+1], c[l+1], a[n+l+1], *dcmplx_a, b[m+l+1], *dcmplx_b, ra[n+l], rb[m+l], rc[l];
  int i, j, c_num, bug, dk, dl, dgcd, k1, k2, dc_gcd;
  dcmplx **gcd_coeff;
  dcmplx *t1, *t2, *c_gcd;  
  int dbd, dad, mult[l];

  printf("please input the polynomial aa\n"); 
  for(i=0; i<=n; i++)
     read_dcmplx( &(aa[i]));
  printf("please input the polynomial bb\n");
  for(i=0; i<=m; i++)
     read_dcmplx(&(bb[i]));
  printf("please input the polynomial c, (the leading coefficient must be one!)\n");
  for(i=0; i<=l; i++)
     read_dcmplx(&(c[i]));

  for(i=0; i<=n+l; i++)
      a[i]=create1(0.0);
  for(i=0; i<=n; i++)
    for(j=0; j<=l; j++)
       a[i+j]=add_dcmplx(a[i+j], mul_dcmplx(aa[i], c[j]));     
  
  for(i=0; i<=m+l; i++)  
    b[i]=create1(0.0);
  for(i=0; i<=m; i++)
    for(j=0; j<=l; j++)
       b[i+j]=add_dcmplx(b[i+j], mul_dcmplx(bb[i], c[j]));

  rootsGCD(n+l+1,a,m+l+1,b,&c_num,ra,rb);

  printf("the roots number of the given common diviser is: %d \n", l);
  multiple_roots(l+1, c, eps, 10*l, rc, mul_tol, mult);
  printf("the root of the gcd is:\n");
  for(i=0; i<l; i++)
     writeln_dcmplx(rc[i]);      
  
  printf("the roots number of the calculated common diviser is: %d \n", c_num); 
  
  if(l!= c_num)
  {
    printf("Bug!!\n");
    exit(1);
  }

  for(i=0; i<l; i++)
  {
    bug=1;
    for(j=0; j<l; j++)
    {
      if(equal_dcmplx(rc[i], ra[j], tol))
        bug=0;
    }
    if( bug==1 )
    {
      printf(" Bug!!\n");
      exit(1);
    }
  }
 printf("the gcd is correct !\n");

 /* now begin to test extended gcd method */
 gcd_coeff=ExtPolyGcd( n+l+1, a, m+l+1, b, &dgcd, &dk, &dl, &dbd, &dad);

 Print_Poly(dl, gcd_coeff[1]);
 printf("a is:\n"); 
 Print_Poly(n+l, a);
 dcmplx_a=(dcmplx*) calloc(n+l+1, sizeof(dcmplx));
 /* for(i=0; i<=n+l; i++)
      dcmplx_a[i]=a[i]; */
  t1=mul_poly( dk, gcd_coeff[0], n+l, a, &k1 );
  t2=mul_poly( dl, gcd_coeff[1], m+l, b, &k2 );


 /*
 printf("t1 is:\n");  
 for(i=0; i<=k1; i++)
   writeln_dcmplx(t1[i]);

 printf("t2 is:\n");
 printf("dl=%d\n", dl);
printf("m+l=%d\n", m+l);
 for(i=0; i<=k2; i++)
    writeln_dcmplx(t2[i]);
 */

 printf("\n"); 
 printf("k is:\n");
 for(i=0; i<=dk; i++)
    writeln_dcmplx(gcd_coeff[0][i]);

 printf("l is:\n");
 for(i=0; i<=dl; i++)
    writeln_dcmplx(gcd_coeff[1][i]);
 printf("\n");

 c_gcd=add_poly( k1, t1, k2, t2, &dc_gcd);

 printf("the given polynomial a is:\n"); 
 for(i=0; i<=n+l; i++)
    writeln_dcmplx(a[i]);
 printf("a's roots are:\n");
 for(i=0; i<n+l; i++)
   writeln_dcmplx(ra[i]);
 printf("\n");

 printf("the given polynomial b is:\n");
 for(i=0; i<=m+l; i++)
    writeln_dcmplx(b[i]);
 printf("b's roots are:\n");
 for(i=0; i<m+l; i++)
   writeln_dcmplx(rb[i]);
 printf("\n");

 printf("the given gcd polynomial is:\n");
 for(i=0; i<=l; i++)
    writeln_dcmplx(c[i]);
 printf("\n");

  if(l!= dc_gcd)
  {
    printf("Bug!! the given gcd is different with the gcd method result.(0 is special case).\n");
    
  }
 
 printf("the calculated gcd polynomial is:\n");
 for(i=0; i<=dc_gcd; i++)
    writeln_dcmplx(c_gcd[i]);
 printf("\n");

 printf("the gcd from common roots method is:\n");
 for(i=0; i<=dgcd; i++)
   writeln_dcmplx(gcd_coeff[2][i]);
 printf("\n");

 free(t1); free(t2);
 free(c_gcd);

 free(gcd_coeff[0]);
 free(gcd_coeff[1]);
 free(gcd_coeff[2]);
 free(gcd_coeff[3]);
 free(gcd_coeff[4]);
 free(gcd_coeff);
}

void random_test ( int n, int m, int l )
{
  double eps = 1.0e-12;
  double tol = 1.0e-8;
  double mul_tol = 1.0e-3;
  dcmplx aa[n+1], bb[m+1], c[l+1], a[n+l+1], b[m+l+1], ra[n+l], rb[m+l], rc[l];
  int i, j, c_num, bug, dk, dl, dgcd, k1, k2, dc_gcd;
  dcmplx **gcd_coeff;
  dcmplx *t1, *t2, *c_gcd;  
  dcmplx near;
  double error=1.0;
  int dbd, dad, mult[l];

  srand(time(NULL));

  for(i=0; i<=n; i++)
    aa[i] = random_dcmplx1();
  for(i=0; i<=m; i++)
    bb[i] = random_dcmplx1();

  c[l]=create1(1.0);
  for(i=0; i<l; i++)
    c[i] = random_dcmplx1();

  for(i=0; i<=n+l; i++)
      a[i]=create1(0.0);
  for(i=0; i<=n; i++)
    for(j=0; j<=l; j++)
       a[i+j]=add_dcmplx(a[i+j], mul_dcmplx(aa[i], c[j]));     
  
  for(i=0; i<=m+l; i++)  
    b[i]=create1(0.0);
  for(i=0; i<=m; i++)
    for(j=0; j<=l; j++)
       b[i+j]=add_dcmplx(b[i+j], mul_dcmplx(bb[i], c[j]));

  rootsGCD(n+l+1,a,m+l+1,b,&c_num,ra,rb);

  printf("the roots number of the given common diviser is: %d \n", l);
  
  multiple_roots(l+1, c, eps, 10*l, rc, mul_tol, mult);
  /*
  printf(" the root of the gcd is:\n");
  for(i=0; i<l; i++)
     writeln_dcmplx(rc[i]);
  printf("\n");      
  */
  printf("the roots number of the calculated common diviser is: %d \n", c_num); 
  
  if(l!= c_num)
  {
    printf("Bug!!\n");
        printf("a is:\n");  
    for(i=0; i<=n+l; i++)
      writeln_dcmplx(a[i]);    
    
    printf("b is:\n");  
    for(i=0; i<=m+l; i++)
      writeln_dcmplx(b[i]);

    printf("a's roots are:\n");
    for(i=0; i<n+l; i++)
      writeln_dcmplx(ra[i]);

    printf("b's roots are:\n");
    for(i=0; i<m+l; i++)
      writeln_dcmplx(rb[i]);

    printf("fail reason: gcd roots number is not correct.\n");
    exit(1);
  }

  for(i=0; i<l; i++)
  {
    bug=1;
    for(j=0; j<l; j++)
    {
      if(equal_dcmplx(rc[i], ra[j], tol))
        bug=0;
      if(dcabs(sub_dcmplx(rc[i], ra[j]))<error)
      {
         error=dcabs(sub_dcmplx(rc[i], ra[j]));
         near=ra[j];        
      }
    }
    if( bug==1 )
    {
     printf("a is:\n");  
     for(i=0; i<=n+l; i++)
      writeln_dcmplx(a[i]);    
    
     printf("b is:\n");  
     for(i=0; i<=m+l; i++)
       writeln_dcmplx(b[i]);

     printf("a's roots are:\n");
     for(i=0; i<n+l; i++)
       writeln_dcmplx(ra[i]);

     printf("b's roots are:\n");
     for(i=0; i<m+l; i++)
       writeln_dcmplx(rb[i]);

      printf("the given common root is:");
      writeln_dcmplx(rc[i]);
      printf("the nearest one to the given common roots is:");
      writeln_dcmplx(near);
      printf("the difference is:%lf\n", error);
   
      printf("fail reason: the common root is not the given gcd root\n");  
      printf(" Bug!!\n");
      exit(1);
    }
  }
 printf("the result is correct !\n");
 
 /* now begin to test extended gcd method */
 gcd_coeff=ExtPolyGcd( n+l+1, a, m+l+1, b, &dgcd, &dk, &dl, &dbd, &dad);
 t1=mul_poly( dk, gcd_coeff[0], n+l, a, &k1 );
 t2=mul_poly( dl, gcd_coeff[1], m+l, b, &k2 );
 
 c_gcd=add_poly( k1, t1, k2, t2, &dc_gcd);

 /*
 printf("the given polynomial a is:\n"); 
 for(i=0; i<=n+l; i++)
    writeln_dcmplx(a[i]);
 printf("\n");

 printf("the given polynomial b is:\n");
 for(i=0; i<=m+l; i++)
    writeln_dcmplx(b[i]);
 printf("\n");
 */

 printf("the given gcd polynomial is:\n");
 for(i=0; i<=l; i++)
    writeln_dcmplx(c[i]);
 printf("\n");
 
 printf("the calculated gcd polynomial is:\n");
 for(i=0; i<=dc_gcd; i++)
    writeln_dcmplx(c_gcd[i]);
 printf("\n");

 printf("the gcd got from common roots is:\n");
 for(i=0; i<=dgcd; i++)
    writeln_dcmplx(gcd_coeff[2][i]);
 printf("\n");

  if(l!= dc_gcd)
  {
    printf("Bug!!\n");

    printf("a is:\n");  
    for(i=0; i<=n+l; i++)
      writeln_dcmplx(a[i]);    
    
    printf("b is:\n");  
    for(i=0; i<=m+l; i++)
      writeln_dcmplx(b[i]);

    printf("a's roots are:\n");
    for(i=0; i<n+l; i++)
      writeln_dcmplx(ra[i]);

    printf("b's roots are:\n");
    for(i=0; i<m+l; i++)
      writeln_dcmplx(rb[i]);

    exit(1);
  }

  if(l!=dgcd)
  {
    printf("Bug!! (get_poly)\n");
    exit(1);
  }

  for(i=0; i<=l; i++)
  {
    bug=1;
    for(j=0; j<=l; j++)
    {
      if(equal_dcmplx(c[i], c_gcd[j], tol))
        bug=0;
    }
    if( bug==1 )
    {
      printf(" Bug!!\n");
      exit(1);
    }
  }

  for(i=0; i<=l; i++)
  {
    bug=1;
    for(j=0; j<=l; j++)
    {
      if(equal_dcmplx(c[i], gcd_coeff[2][j], tol))
        bug=0;
    }
    if( bug==1 )
    {
      printf(" Bug!!(get_poly)\n");
      exit(1);
    }
  }
 printf(" the extended gcd result is correct !\n");

 free(t1); free(t2);
 free(c_gcd);

 free(gcd_coeff[0]);
 free(gcd_coeff[1]);
 free(gcd_coeff[2]);
 free(gcd_coeff);
}


