#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poly_matrix.h"
#include "poly_gcd.h"
#include "poly_smith.h"

#define  tol 10e-8

void Find_pivots (int n, int m, POLY a[n][m], int * pivrow, int *pivcol)
{
  int i,j, min, row, col;
  int  min_degree=20; 

  row=*pivrow;
  col=*pivcol;

  for ( i=*pivrow; i<n; i++)
      for ( j=*pivcol; j<m; j++ )
  
	if(((a[i][j].d>0) || (!equal_dcmplx( a[i][j].p[0], zero, tol))) && (a[i][j].d<min_degree)) 
	{
	  min_degree=a[i][j].d; 
	  row=i;
	  col=j;
	}
   *pivrow=row;
   *pivcol=col;

}

void Permute ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m],
			   int pivrow, int row, int pivcol, int col)
{
    if (pivrow != row)
    {
      Interchange_Rows ( n, m, a, pivrow, row);
      Interchange_Rows ( n, n, p, pivrow, row);
    }
    if (pivcol != col)
    {
      Interchange_Cols ( n, m, a, pivcol, col);
      Interchange_Cols ( m, m, q, pivcol, col);
    }
}

void Interchange_Rows ( int n, int m, POLY a[n][m], int r1, int r2)
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

void Interchange_Cols ( int n, int m, POLY a[n][m], int c1, int c2)
{
  int i;
  POLY temp;
  for ( i=0; i<n; i++)
    {
       temp=a[i][c1];
       a[i][c1]=a[i][c2];
       a[i][c2]=temp;
    }
}


void Eliminate_Col ( int n, int m, POLY a[n][m], POLY p[n][n], int r, int c )
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

void Eliminate_Row ( int n, int m, POLY a[n][m], POLY q[m][m], int r, int c )
{
  int i, j, dgcd, dk, dl, dbd, dad;
  POLY   aa, bb, neg_bd, tmp, t1, t2, bd, ad;
  dcmplx **gcd_coeff;
  dcmplx *d;
 
  for (j=c+1; j<m; j++)
  {
     aa=assign_poly(a[r][c]);
     bb=assign_poly(a[r][j]);

     if (!equal_dcmplx( bb.p[0], zero, tol) || (bb.d>0))
     {
        gcd_coeff=ExtPolyGcd( aa.d+1, aa.p, bb.d+1, bb.p, &dgcd, &dk, &dl, &dbd, &dad);   
        d = gcd_coeff[2];
        bd.p = gcd_coeff[3]; bd.d=dbd;
        neg_bd = neg_poly(bd);
        ad.p = gcd_coeff[4]; ad.d=dad;

        for (i=0; i<n; i++)
        {  
           tmp = assign_poly(a[i][c]);
           t1.p = mul_poly( dk, gcd_coeff[0], tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( dl, gcd_coeff[1], a[i][j].d, a[i][j].p, &(t2.d));
           free(a[i][c].p);        
       
           a[i][c].p = add_poly( t1.d, t1.p, t2.d, t2.p, &(a[i][c].d));
           free(t1.p); free(t2.p);
           /* a[i][c]=x*tmp+y*a[i][j];*/

           t1.p = mul_poly( neg_bd.d, neg_bd.p, tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( ad.d, ad.p, a[i][j].d, a[i][j].p, &(t2.d));
           free(a[i][j].p);
           a[i][j].p = add_poly(t1.d, t1.p, t2.d, t2.p, &(a[i][j].d));
           free(t1.p); free(t2.p);
           free(tmp.p);
           /* a[i][j]=(-bb/d)*tmp+(aa/d)*a[i][j]; */

        }
 		   
        for (i=0; i<m; i++)
        { 
           tmp = assign_poly(q[i][c]);

           t1.p = mul_poly( dk, gcd_coeff[0], tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( dl, gcd_coeff[1], q[i][j].d, q[i][j].p, &(t2.d));
           free(q[i][c].p);

           q[i][c].p = add_poly( t1.d, t1.p, t2.d, t2.p, &(q[i][c].d));
           free(t1.p); free(t2.p);   
           /* q[i][c]=x*tmp+y*q[i][j]; */

           t1.p = mul_poly( neg_bd.d, neg_bd.p, tmp.d, tmp.p, &(t1.d));
           t2.p = mul_poly( ad.d, ad.p, q[i][j].d, q[i][j].p, &(t2.d));
           free(q[i][j].p);
           q[i][j].p = add_poly(t1.d, t1.p, t2.d, t2.p, &(q[i][j].d));
           free(t1.p); free(t2.p);
           free(tmp.p);
           /* q[i][j]=(-bb/d)*tmp+(aa/d)*q[i][j]; */
        }
	free(gcd_coeff[0]); free(gcd_coeff[1]); free(gcd_coeff[2]); free(gcd_coeff);
        free(bd.p); free(ad.p); free(neg_bd.p);
     }
     free(aa.p); free(bb.p);
  }	

}


void Smith_Diagonal ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m])
{
    int row, col, pivrow, pivcol;

    while (Diagonal ( n, m , a) !=1 )
	{ 
		row=0;
		col=0;
	  while ((row<=n-1) && (col<=m-1))   /* different with before (&&) */
	  {
	    pivrow=row;
	    pivcol=col;

  	    Find_pivots ( n, m, a, &pivrow, &pivcol);
            
	    if( (a[pivrow][pivcol].d==0) && equal_dcmplx(a[pivrow][pivcol].p[0],zero,tol)) 
	    {
               break;
            }

            Permute ( n, m, a, p, q, pivrow, row, pivcol, col); 
	    Eliminate_Col ( n, m, a, p, row, col );
	    Eliminate_Row ( n, m, a, q, row, col );

	    row++;
	    col++;
	  }

	}


}

int Diagonal ( int n, int m , POLY a[n][m])
{
	int i, j;
	for ( i=0; i<n; i++)
		for ( j=0; j<m; j++)
			if ( i!=j )
			  if((a[i][j].d>0) || (!equal_dcmplx( a[i][j].p[0], zero, tol)))
					return 0;

    return 1;
}    

int poly_divide ( POLY a, POLY b )
{
  int i, dgcd, dk, dl, dbd, dad;
  dcmplx **gcd_coeff;

  gcd_coeff=ExtPolyGcd( a.d+1, a.p, b.d+1, b.p, &dgcd, &dk, &dl, &dbd, &dad);
  if(((dl == 0) && equal_dcmplx( gcd_coeff[1][0], zero, tol)) ||
     ((b.d== 0) && equal_dcmplx(b.p[0], zero, tol)))
  { 

    free_gcd_coeff(gcd_coeff);
    return 1; 
  }

  free_gcd_coeff(gcd_coeff);
  return 0;

}

void add_row2above(int n, int m, POLY a[n][m], POLY p[n][n], int row)
{
  int j;
  POLY tmp;

  a[row-1][row]=assign_poly(a[row][row]);    
  for ( j=0; j<n; j++ )
  { 
    tmp =add_poly1 (p[row-1][j], p[row][j]);
    free(p[row-1][j].p);
    p[row-1][j]=assign_poly(tmp); 
   /* free(tmp); */

  }

}

void Smith ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m])
{
  int i, min;
  
  I_matrix ( n, p );
  I_matrix ( m, q );
  Smith_Diagonal( n, m, a, p, q );

  if(n<=m) min=n;
  else min=m;

  i = 0;
  while(i<min-1)
  {  
    if(!poly_divide(a[i][i], a[i+1][i+1]))
    { 
       add_row2above(n, m, a, p, i+1);
       Smith_Diagonal( n, m, a, p, q );
       i=0;
    }   

    i++;
  } 
}






