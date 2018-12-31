#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "dcmplx.h"
#include "dc_matrix.h"
#include "poly_matrix.h"
#include "poly_gcd.h"
#include "poly_hermite.h"

#define  tol 10e-8

int pivot_row (int n, int m, POLY a[n][m], int k, int z)
{
  int i, j, min, row, flag;
  int  min_degree=20; 

  flag=0;
  j=k+z;

  for ( i=k; i<n; i++) 
	if(((a[i][j].d>0) || (!equal_dcmplx( a[i][j].p[0], zero, tol))) && (a[i][j].d<min_degree)) 
	{ 
          flag = 1;
	  min_degree = a[i][j].d; 
	  row = i;
	}
  if(flag==1)
    return row;
  else
    return -1;   

}

void Permute_row ( int n, int m, POLY a[n][m], POLY p[n][n], int pivrow, int row)
{
    if (pivrow != row)
    {
      Interchange_Rows1 ( n, m, a, pivrow, row);
      Interchange_Rows1 ( n, n, p, pivrow, row);
    }
}

void Interchange_Rows1 ( int n, int m, POLY a[n][m], int r1, int r2)
{
  int j;
  POLY temp;

  for ( j=0; j<m; j++)
    {
       temp=a[r1][j];
       a[r1][j]=a[r2][j];
       a[r2][j]=temp;
    }
}

void Eliminate_Col1 ( int n, int m, POLY a[n][m], POLY p[n][n], int r, int c )
{
  int i, j, dgcd, dk, dl, dbd, dad;
  POLY   aa, bb, neg_bd, tmp, t1, t2, bd, ad;
  dcmplx **gcd_coeff;
  dcmplx *d;

  for (i=r+1; i<n; i++)
  {
     aa=assign_poly(a[r][c]);
     bb=assign_poly(a[i][c]);

     if (!equal_dcmplx( bb.p[0], zero, tol) || (bb.d>0))
     { 
        gcd_coeff=ExtPolyGcd( aa.d+1, aa.p, bb.d+1, bb.p, &dgcd, &dk, &dl, &dbd, &dad);   
	/*printf("aa="); Print_Poly(aa.d, aa.p);
	printf("bb=");Print_Poly(bb.d,bb.p);
        printf("gcd=");Print_Poly(dgcd, gcd_coeff[2]);
        */
 
        d = gcd_coeff[2];
        bd.p = gcd_coeff[3]; bd.d=dbd;
        neg_bd = neg_poly(bd);
        ad.p = gcd_coeff[4]; ad.d=dad;
       
        for (j=0; j<m; j++)
        { 
           tmp = assign_poly(a[r][j]);
	   t1.p = mul_poly( dk, gcd_coeff[0], tmp.d, tmp.p, &(t1.d));
	   t2.p = mul_poly( dl, gcd_coeff[1], a[i][j].d, a[i][j].p, &(t2.d));
           free(a[r][j].p);
           a[r][j].p = add_poly( t1.d, t1.p, t2.d, t2.p, &(a[r][j].d));
           free(t1.p); free(t2.p);
	   /* a[r][j]=x*tmp+y*a[i][j];*/
	   /* printf("a[%d][%d]=",r, j); Print_Poly(a[r][j].d, a[r][j].p); */          
           t1.p = mul_poly( neg_bd.d, neg_bd.p, tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( ad.d, ad.p, a[i][j].d, a[i][j].p, &(t2.d));
           free(a[i][j].p);      
           a[i][j].p = add_poly(t1.d, t1.p, t2.d, t2.p, &(a[i][j].d));       

           free(t1.p); free(t2.p);
           free(tmp.p);
        /* a[i][j]=(-bb/d)*tmp+(aa/d)*a[i][j];*/
        }
	
	 
        for (j=0; j<n; j++)
        {
           tmp =assign_poly(p[r][j]);
           t1.p = mul_poly( dk, gcd_coeff[0], tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( dl, gcd_coeff[1], p[i][j].d, p[i][j].p, &(t2.d));
           free(p[r][j].p);

           p[r][j].p = add_poly( t1.d, t1.p, t2.d, t2.p, &(p[r][j].d));
           free(t1.p); free(t2.p);      
           /* p[r][j]=x*tmp+y*p[i][j]; */

           t1.p = mul_poly(neg_bd.d, neg_bd.p, tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( ad.d, ad.p, p[i][j].d, p[i][j].p, &(t2.d));
           free(p[i][j].p);
           p[i][j].p = add_poly(t1.d, t1.p, t2.d, t2.p, &(p[i][j].d));
           free(t1.p); free(t2.p);
           free(tmp.p);
           /* p[i][j]=(-bb/d)*tmp+(aa/d)*p[i][j]; */

        }
	free(gcd_coeff[0]); free(gcd_coeff[1]); free(gcd_coeff[2]); free(gcd_coeff);
        free(bd.p); free(ad.p);  free(neg_bd.p);
     }
     free(aa.p); free(bb.p);
  }	
}

void Hermite ( int n, int m, POLY a[n][m], POLY p[n][n])
{
  int i, j, k, z;
  
  I_matrix ( n, p );
  z=0;  

  for(j=0; j<m; j++)
  {
    k=j-z;
    if(k<n-1)
    {
      i=pivot_row(n, m, a, k, z);

      if(i==-1)
        z++;
      else
      {
        if(i!=k) 
	  Permute_row ( n, m, a, p, k, i);
	Eliminate_Col1( n, m, a, p, k, j);
      }
    }
    /* printf("Here is the matrix in %d the iterative.\n", j); print(n, m, a); */
  }
}






