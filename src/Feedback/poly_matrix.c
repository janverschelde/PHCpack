#include <stdlib.h>
#include "poly_matrix.h"
#include "poly_gcd.h"
#include "poly_smith.h"
#include "dc_roots.h"

void read_matrix ( int n, int m, POLY a[n][m] )
{
   int i, j, k;
   printf("Give entries of %d-by-%d matrix :\n", n, m);
   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
      {
         printf("Please input the degree of the Polynomial ");
         printf("( %d, %d )\n", i, j);
         scanf("%d",&(a[i][j].d));
         a[i][j].p = (dcmplx*) calloc(a[i][j].d+1, sizeof(dcmplx)); 
         printf("Please input the polynomial:\n");
         for (k=0; k<=a[i][j].d; k++) 
	   read_dcmplx(&(a[i][j].p[k]));
         
      }
}

void read_matrix1 ( int n, int m, POLY a[n][m] )
{
   int i, j, k;
   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
      {
         scanf("%d",&(a[i][j].d));
         a[i][j].p = (dcmplx*) calloc(a[i][j].d+1, sizeof(dcmplx));
         for (k=0; k<=a[i][j].d; k++)
           read_dcmplx(&(a[i][j].p[k]));

      }
}

void random_matrix ( int n, int m, POLY a[n][m] )
{
	int i, j;
        POLY gcd[n][m];
        POLY tmp;   
        
        for(i=0; i<n; i++)  
          for(j=0; j<m; j++)
            gcd[i][j]=random_poly1(i);

        tmp=random_poly1(1);
	for (i=0; i<n; i++)
		for (j=0; j<m; j++)
		  {
		     a[i][j]=mul_poly1(tmp, gcd[i][j]);                     
		  }

        for(i=0; i<n; i++)  
          for(j=0; j<m; j++)
            free(gcd[i][j].p);

        free(tmp.p); 

     
}

void I_matrix ( int n, POLY p[n][n])
{
  int i, j;
 
  for ( i=0; i<n; i++)
     for ( j=0; j<n; j++)
        { 
          if ( i==j )
	  { 
            p[i][j].p = (dcmplx*) calloc(1, sizeof(dcmplx));
            p[i][j].d = 0;
            p[i][j].p[0]=one;
          }
          else 
	  {           
            p[i][j].p = (dcmplx*) calloc(1, sizeof(dcmplx));
            p[i][j].d = 0;
            p[i][j].p[0]=zero;
          }  
        }
}


void zero_matrix ( int n, int m, POLY p[n][m])
{
  int i, j;

  for ( i=0; i<n; i++)
     for ( j=0; j<m; j++)
        { 
          p[i][j].p = (dcmplx*) calloc(1, sizeof(dcmplx));
          p[i][j].d = 0;
          p[i][j].p[0]=zero;  
        }
}

void copy ( int n, int m, POLY a[n][m], POLY t[n][m])
{
   int i,j;
   for (i=0; i<n; i++)
       for (j=0; j<m; j++)
           t[i][j]=assign_poly(a[i][j]);
 
}

void print ( int n, int m, POLY a[n][m] )
{
   int i,j;
   for (i=0; i<n; i++)
   {
      for (j=0; j<m; j++)
      {
        printf("( %d, %d )\n", i, j); 
        Print_Poly(a[i][j].d, a[i][j].p);
        printf("\n");
      }
   } 
}

void print1 ( int n, int m, POLY a[n][m] )
{
   int i,j;
   for (i=0; i<n; i++)
   {
      for (j=0; j<m; j++)
      {
	/* printf("( %d, %d )\n", i, j); */
        printf("%d\n", a[i][j].d);
        Print_Poly(a[i][j].d, a[i][j].p);
        printf("\n");
      }
   }
}

void Multiply ( int n, int m, int l, POLY a[n][m], POLY b[m][l], POLY c[n][l] )
{
    int i, j, k;
    POLY tmp1, tmp2;

    zero_matrix (  n, l, c);


    for ( i=0; i<n; i++ )
	{
      for ( j=0; j<l; j++ )
	  {
          for ( k=0; k<m; k++ )
		  {

                     tmp1.p=mul_poly(a[i][k].d, a[i][k].p, b[k][j].d, b[k][j].p, &(tmp1.d));

                     tmp2.p=add_poly(c[i][j].d, c[i][j].p, tmp1.d, tmp1.p, &(tmp2.d));
	             free(c[i][j].p);

                     c[i][j]=assign_poly(tmp2);
                     free(tmp1.p); 
		     free(tmp2.p);
		    
		    /* c[i][j]+=a[i][k]*b[k][j];*/

		  }
	  }
	}
}


void Transpose ( int n, int m, POLY a[n][m], POLY b[m][n] )
{
	int i,j;
	for ( i=0; i<m; i++ )
		for ( j=0; j<n; j++)
			b[i][j]=assign_poly(a[j][i]);
}


void free_matrix ( int n, int m, POLY a[n][m])
{
     int i, j;
     for ( i=0; i<n; i++ )
          for ( j=0; j<m; j++)
              free(a[i][j].p);
}

POLY Inverse_Poly ( int n, POLY (*M)[n] )
{ 
   POLY P[n][n], Q[n][n], t_M1[n][n], t_M2[n][n];
   POLY ds;
   int i, j, dk, dl, dgcd, dc_gcd, dad, dbd;
   dcmplx **gcd_coeff, tmp;

   copy(n, n, M, t_M1);  

   Smith(n, n, t_M1, P, Q);
   /* printf("the smith form is:\n");  print1(n, n, t_M1); */
   ds=assign_poly(t_M1[n-1][n-1]);

   for( i=0; i<n-1; i++ )
     { 
       gcd_coeff = ExtPolyGcd(ds.d+1, ds.p, t_M1[i][i].d+1, t_M1[i][i].p, &dgcd, &dk, &dl, &dbd, &dad);  
       free(t_M1[i][i].p);
       t_M1[i][i].d = dad;
       t_M1[i][i].p = assign(dad,gcd_coeff[4]);
       for( j=0; j<5; j++)
         free(gcd_coeff[j]);
       free(gcd_coeff);
     }
   free(t_M1[n-1][n-1].p);

   t_M1[n-1][n-1].p = (dcmplx*) calloc(1, sizeof(dcmplx));
   t_M1[n-1][n-1].d = 0;
   t_M1[n-1][n-1].p[0] = one;

   /* calculate ds*Q*(t_M1's inverse) * P = ds*M  */
   Multiply( n, n, n, Q, t_M1, t_M2 );

   free_matrix ( n, n, M );
   Multiply( n, n, n, t_M2, P, M );

   if(ds.p[ds.d].re<0)
   {
     negative(ds.d, ds.p);
     neg_polymatrix( n, n, M );
   }   

   /* make the leading coefficient of ds one */
   tmp = ds.p[ds.d];
   divide_by_number(ds.d, ds.p, tmp); 

   for(i=0; i<n; i++)
     for(j=0; j<n; j++)
       {  
        divide_by_number(M[i][j].d, M[i][j].p, tmp);
       }
   /* printf("ds="); Print_Poly(ds.d , ds.p); */
  
   free_matrix ( n, n, t_M1 );
   free_matrix ( n, n, t_M2 );
   free_matrix ( n, n, P );
   free_matrix ( n, n, Q );
   return ds;
}


void dcmatrix_Multiply ( int n, int m, int l, dcmplx a[n][m], POLY b[m][l], POLY c[n][l] )
{
    int i, j, k;
    POLY tmp1, tmp2;

    zero_matrix (  n, l, c);
    for ( i=0; i<n; i++ )
	{
      for ( j=0; j<l; j++ )
	  {
          for ( k=0; k<m; k++ )
		  {

                     tmp1 = mul_dcmplx_poly ( a[i][k], b[k][j] );
                     tmp2.p=add_poly(c[i][j].d, c[i][j].p, tmp1.d, tmp1.p, &(tmp2.d));
	             free(c[i][j].p);
                     c[i][j]=assign_poly(tmp2);
                     free(tmp1.p); 
		     free(tmp2.p);
		    
		    /* c[i][j]+=a[i][k]*b[k][j];*/

		  }
	  }
	}

}

void add_polymatrix( int n, int m, POLY a[n][m], POLY b[n][m], POLY c[n][m] )
{
  int i, j;
  for(i=0; i< n; i++ )
    for(j=0; j<m; j++ )
       c[i][j] = add_poly1(a[i][j], b[i][j]); 

}

void sub_polymatrix( int n, int m, POLY a[n][m], POLY b[n][m], POLY c[n][m] )
{
  int i, j;
  for(i=0; i< n; i++ )
    for(j=0; j<m; j++ )
       c[i][j] = min_poly1(a[i][j], b[i][j]); 

}


void evaluate_matrix( int n, int m, POLY a[n][m], dcmplx b[n][m], dcmplx x)
{
  int i, j;
  
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
    { 
      if(a[i][j].d>=1)
        b[i][j] = horner(a[i][j].d+1, a[i][j].p, x);
      else
        b[i][j] = a[i][j].p[0];
    }
}



void neg_polymatrix( int n, int m, POLY a[n][m] )
{
  int i, j;
  
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
      negative( a[i][j].d, a[i][j].p );
}
 
