#include <stdio.h>
#include <stdlib.h>
#include "poly_dcmplx.h"
#include "dcmplx.h"
#define tol 10e-8

void Read_Poly1 ( int n, dcmplx a[] )
{
  int i;
  for (i=0; i<=n; i++)
	  read_dcmplx(&(a[i]));
}
 

void Read_Poly ( int n, dcmplx a[], int m, dcmplx b[] )
{
  int i, j;
  printf("Please input the first polynomial:\n");
  for (i=0; i<=n; i++)
	  read_dcmplx(&(a[i]));
  printf("Please input the second polynomial:\n");
  for (j=0; j<=m; j++)
	  read_dcmplx(&(b[j]));
  
  printf("The polynomial a is:\n");
  Print_Poly(n,a);
  printf("The polynomial b is:\n");
  Print_Poly(m,b);
}

POLY random_poly (int max_deg)
{
  POLY a;
  int i;

  a.d = (int) rand()%(max_deg+1);
  a.p = (dcmplx*) calloc(a.d+1, sizeof(dcmplx));
  for ( i=0; i<=a.d; i++ )
      a.p[i] = random_dcmplx1();
  return a;
}  

POLY random_poly1 (int deg)
{
  POLY a;
  int i;

  a.d =  deg;
  a.p = (dcmplx*) calloc(a.d+1, sizeof(dcmplx));
  for ( i=0; i<=a.d; i++ )
      a.p[i] = random_dcmplx1();
  return a;
}  



void Print_Poly ( int k, dcmplx * c )
{
	int i;
	for ( i=0; i<=k; i++ )
	{
	 writeln_dcmplx(c[i]);
	}
}

POLY neg_poly ( POLY a)
{
  int i;
  POLY b;
  b.p = (dcmplx*) calloc(a.d+1, sizeof(dcmplx));
  for ( i=0; i<=a.d; i++)
      b.p[i] = min_dcmplx(a.p[i]);
  b.d = a.d;
  return b;
}


dcmplx* add_poly ( int n, dcmplx a[], int m, dcmplx b[], int *k )
{
  int i;
  dcmplx *c;

  if ( m>=n ) 
  {
      *k=m;
      c = (dcmplx*) calloc(*k+1, sizeof(dcmplx));

      for ( i=0; i<=n; i++ )
	    c[i]=add_dcmplx(a[i], b[i]);
      for ( i=n+1; i<=m; i++ )
	    c[i]=b[i]; 

  }
  else 
  {
      *k=n;
      c = (dcmplx*) calloc(*k+1, sizeof(dcmplx));
 
      for ( i=0; i<=m; i++ )
	    c[i]=add_dcmplx(a[i], b[i]);
      for ( i=m+1; i<=n; i++ )
	     c[i]=a[i];
  }

  *k=degree( c, *k );  

  return c;
}

POLY add_poly1 ( POLY a, POLY b)
{
  POLY c;
  c.p = add_poly(a.d, a.p, b.d, b.p, &(c.d));
  return c;
}


dcmplx* min_poly ( int n, dcmplx a[], int m, dcmplx b[], int *k )
{
  int i;
  dcmplx *c;

  if ( m>=n ) 
  {
      *k=m;
      c = (dcmplx*) calloc(*k+1, sizeof(dcmplx));
      for ( i=0; i<=n; i++ )
	    c[i]=sub_dcmplx(a[i], b[i]);
      for ( i=n+1; i<=m; i++ )
	    c[i]=sub_dcmplx(c[i], b[i]); 
  }
  else 
  {
      *k=n;
      c = (dcmplx*) calloc(*k+1, sizeof(dcmplx));
 
      for ( i=0; i<=m; i++ )
	    c[i]=sub_dcmplx(a[i], b[i]);
      for ( i=m+1; i<=n; i++ )
	     c[i]=a[i];
  }
    *k=degree( c, *k ); 
  return c;
}

POLY min_poly1( POLY a, POLY b)
{
  POLY c;
  c.p = min_poly( a.d, a.p, b.d, b.p, &(c.d));
  return c;
}
dcmplx* mul_poly ( int n, dcmplx a[], int m, dcmplx b[], int *k )
{
  int i, j;
  dcmplx *c; 
  if (iszero( n,a ) || iszero ( m, b ) )
  {
    *k=0;
    c = (dcmplx*) calloc(1, sizeof(dcmplx));
    c[0]=create1(0.0);
    return c;
  }
  *k=n+m;
  c = (dcmplx*) calloc(*k+1, sizeof(dcmplx));

  for ( i=0; i<=n; i++)
	  for ( j=0; j<=m; j++)
	    {
	      c[i+j]=add_dcmplx(c[i+j],mul_dcmplx(a[i],b[j]));            
            }
  *k=degree( c, *k ); 
  return c;
}

POLY mul_poly1 ( POLY a, POLY b)
{
  POLY c;
  c.p=mul_poly ( a.d, a.p, b.d, b.p, &(c.d));
  return c;
}

POLY mul_dcmplx_poly ( dcmplx a, POLY b )
{
  POLY c;
  int i;
  c.d = b.d;
  c.p = (dcmplx*) calloc(c.d+1, sizeof(dcmplx));
  for(i=0; i<=c.d; i++)
    c.p[i] = mul_dcmplx( a, b.p[i]);
  return c; 

}
dcmplx** div_poly ( int n, dcmplx a[], int m, dcmplx b[],
			    int* dq, int* dr )
{
  int i, j, k, dtemp, d;
  dcmplx *temp_a, *temp_b;
  dcmplx ** result=(dcmplx**) calloc(2, sizeof(dcmplx *));
  if ( iszero (m, b) ) 
  {
	  printf("The divisor can not be 0!\n"); 
	  return result;
  }
  if ( m>n )
  {
	  *dq=0;
	  result[0] = (dcmplx*) calloc(*dq+1, sizeof(dcmplx));
	  result[0][0]=create1(0.0);

	  *dr=n;
	  result[1] = (dcmplx*) calloc(*dr+1, sizeof(dcmplx));
	  for ( i=0; i<=*dr; i++)
		  result[1][i]=a[i];
  }
  else
  {
	  dtemp=n;
	  temp_a = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
	  for ( i=0; i<=n; i++ )
		  temp_a[i]=a[i];
	  *dq=n-m;
	  result[0] = (dcmplx*) calloc(*dq+1, sizeof(dcmplx));
	  while ( dtemp>=m && (!iszero( dtemp,temp_a )))
	  {  
              d=dtemp-m;
	      result[0][d]=div_dcmplx(temp_a[dtemp], b[m]);
	
	      temp_b = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
              for ( j=0; j<=m; j++ )
			  temp_b[d+j]=mul_dcmplx(result[0][d], b[j]);
	      for ( k=0; k<=dtemp; k++ )
			  temp_a[k]=sub_dcmplx(temp_a[k],temp_b[k]);
 	 
	      dtemp= degree ( temp_a, dtemp );           
               
	      free(temp_b);
          }

	  *dr=dtemp;
          result[1] = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
	  for ( i=0; i<=dtemp; i++ )
			result[1][i]=temp_a[i];
          free(temp_a);
	
  } 

        return result;
}

dcmplx* div_poly1 ( int n, dcmplx a[], int m, dcmplx b[], int *dq )
{
  int i, j, k, dtemp, d;
  dcmplx *temp_a, *temp_b;
  dcmplx * result;
  if ( iszero (m, b) ) 
  {
	  printf("The divisor can not be 0!\n"); 
	  return result;
  }
 
  /* n>=m */
   dtemp=n;  printf("dtemp=%d\n", dtemp);
 printf("ok2\n");   temp_a = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
  printf("ok1\n");  for ( i=0; i<=n; i++ )
        temp_a[i]=a[i];
   *dq=n-m;
   result = (dcmplx*) calloc(*dq+1, sizeof(dcmplx));
   while ( dtemp>=m && (!iszero( dtemp,temp_a )))
	  {  
              d=dtemp-m;
	      result[d]=div_dcmplx(temp_a[dtemp], b[m]);
	
	      temp_b = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
              for ( j=0; j<=m; j++ )
			  temp_b[d+j]=mul_dcmplx(result[d], b[j]);
	      for ( k=0; k<=dtemp; k++ )
			  temp_a[k]=sub_dcmplx(temp_a[k],temp_b[k]);
 	 
	      dtemp= degree ( temp_a, dtemp );           
               
	      free(temp_b);
          }

    free(temp_a); 
printf("***\n");
    return result;
}



dcmplx* div_poly2 ( int n, dcmplx a[], int m, dcmplx b[],
			    int* dq )
{
  int i, j, k, dtemp, d;
  dcmplx *temp_a, *temp_b;
  dcmplx ** result=(dcmplx**) calloc(2, sizeof(dcmplx *));
  if ( iszero (m, b) ) 
  {
	  printf("The divisor can not be 0!\n"); 
	  return result[0];
  }
  if ( m>n )
  {
	  *dq=0;
	  result[0] = (dcmplx*) calloc(*dq+1, sizeof(dcmplx));
	  result[0][0]=create1(0.0);

	  /*  *dr=n;
	      result[1] = (dcmplx*) calloc(*dr+1, sizeof(dcmplx));
	      for ( i=0; i<=*dr; i++)
	        result[1][i]=a[i]; */
  }
  else
  {
	  dtemp=n;
	  temp_a = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
	  for ( i=0; i<=n; i++ )
		  temp_a[i]=a[i];
	  *dq=n-m;
	  result[0] = (dcmplx*) calloc(*dq+1, sizeof(dcmplx));
	  while ( dtemp>=m && (!iszero( dtemp,temp_a )))
	  {  
              d=dtemp-m;
	      result[0][d]=div_dcmplx(temp_a[dtemp], b[m]);
	
	      temp_b = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
              for ( j=0; j<=m; j++ )
			  temp_b[d+j]=mul_dcmplx(result[0][d], b[j]);
	      for ( k=0; k<=dtemp; k++ )
			  temp_a[k]=sub_dcmplx(temp_a[k],temp_b[k]);
 	 
	      dtemp= degree ( temp_a, dtemp );           
               
	      free(temp_b);
          }

	  /* *dr=dtemp;
	     result[1] = (dcmplx*) calloc(dtemp+1, sizeof(dcmplx));
	     for ( i=0; i<=dtemp; i++ )
         	     result[1][i]=temp_a[i]; */
	  free(temp_a);
	
  } 

        return result[0];
}

int degree ( dcmplx *a, int d )
{
  int i;
  for ( i=d; i>0; i-- )
      if ( !equal_dcmplx( a[i], create1(0), tol ))
         break;
  return i;
}

int iszero ( int n,  dcmplx a[] )
{
  int i;
  for ( i=0; i<=n; i++ )
      if (! equal_dcmplx(a[i], create1(0), tol) )
      return 0;
  return 1;
}
/*
double dabs (double x)
{
  if (x>=0.0) return x;
  return -x;
}
*/
/*
void Test_Div ( int n, int m)
{
	int i, j, k, num;
	dcmplx a[n+1], b[m+1];

	int dq, dr, kk;
    dcmplx * t_aa, *aa;
    dcmplx ** qr;

    printf("Please input how many tests you want.\n");
	scanf("%d", &num);
	srand(time(NULL));
	for ( k=1; k<=num; k++ )
	{
      printf("\n*******Testing the No. %d case.***********\n", k);
      for ( i=0; i<=n; i++ )
		  a[i]= ((dcmplx) rand())/RAND_MAX;
	
	  for ( j=0; j<=m; j++ )
          b[j]=((dcmplx) rand())/RAND_MAX;
      
      qr=div_poly ( n, a, m, b, &dq, &dr );

	  Print_Poly ( n, a );
	  printf("\t="); Print_Poly ( m, b );
	  printf(" \t* "); Print_Poly ( dq, qr[0] );
	  printf(" \t+ "); Print_Poly ( dr, qr[1] ); 

      t_aa=mul_poly ( dq, qr[0], m, b, &kk ); 
      aa=add_poly ( kk, t_aa, dr, qr[1], &kk );
     

	  free(t_aa);
      free(qr[0]);
      free(qr[1]);
	  free(qr);
	  if ( kk != n ) 
	  {
		
		printf(" BUG!!! (dgree not equal)\n");
		exit(1);
	  } 
	  for ( i=0; i<n; i++ )
		if ( dabs(a[i] - aa[i])>tol )
		{
          printf("i=%d,%12.10lf\n", i,dabs(a[i] - aa[i]));
		  printf(" BUG!!!\n");
		  exit(1);
        }
      printf("OK!\n");
	  free(aa);
 
	}
}
*/
dcmplx *assign(int n, dcmplx a[])
{
    int i;
    dcmplx* b=(dcmplx*) calloc(n+1, sizeof(dcmplx));
    for(i=0; i<=n; i++)
    b[i]=a[i];
    return b;
}

POLY assign_poly(POLY a)
{
    POLY b;
    b.p = assign(a.d, a.p);
    b.d = a.d;
    return b;
}   


void negative(int n, dcmplx a[])
{      int i;
       if(!iszero(n, a))
          for(i=0; i<=n; i++)
              a[i]=min_dcmplx(a[i]);
}        

int equal_poly(int n, dcmplx a[], int m, dcmplx b[])
{
	dcmplx *c;
        int k;
	c= min_poly ( n, a, m, b, &k );
        if ( iszero( k, c ) )
        {  
           free(c); 
           return 1;
         } 
             
        else 
        {
            free(c); 
            return 0;
         }
}

void divide_by_number( int n, dcmplx a[], dcmplx num )
{ 
        int i;
        for(i=0; i<=n; i++)
            a[i]=div_dcmplx(a[i],num);
}

void mult_by_number( int n, dcmplx a[], dcmplx num)
{
        int i; 
        for(i=0; i<=n; i++)
	  a[i]=mul_dcmplx(a[i], num);
}
