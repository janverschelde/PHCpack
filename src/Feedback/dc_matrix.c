#include "dc_matrix.h"
#include <stdlib.h>
#include <math.h>

void print_dcmatrix ( int n, int m, dcmplx a[n][m] )
{
   int i,j;
   for (i=0; i<n; i++)
   {
      for (j=0; j<m; j++)
      {
        printf("( %d, %d )\n", i, j); 
        writeln_dcmplx(a[i][j]);
      }
   } 
}


void print_dcmatrix1 ( int n, int m, dcmplx a[n][m], FILE *ofp)
{
   int i,j;
   for (i=0; i<n; i++)
   {
      for (j=0; j<m; j++)
      {
	// fprintf(ofp,"( %d, %d )\n", i, j); 
        write_dcmplx1(a[i][j], ofp);
        fprintf(ofp, "      ");
      }
      fprintf(ofp, "\n");
   } 
}

void read_dcmatrix ( int n, int m, dcmplx a[n][m] )
{
   int i,j;

   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
      {
         printf("Give value for [%d][%d] : \n", i, j);
         read_dcmplx(&a[i][j]);
      }
}

void read_dcmatrix0 ( int n, int m, dcmplx a[n][m], FILE *ifp )
{
  int i, j;

  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      read_dcmplx0(&a[i][j], ifp);

}

void read_dcmatrix1 ( int n, int m, dcmplx a[n][m] )
{
   int i,j;

   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
      {
         read_dcmplx(&a[i][j]);
      }
}

void read_dcmatrix2 ( int n, int m, dcmplx a[n][m], FILE *ifp )
{
   int i,j;

   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
      {
         read_dcmplx1(&a[i][j], ifp);
      }
}


void random_dcmatrix ( int n, int m, dcmplx a[n][m] )
{
   int i,j;

   for(i=0; i<n; i++)
      for(j=0; j<m; j++)
         a[i][j] = random_dcmplx1();
}

void random_dcmatrix0 ( int n, int m, dcmplx a[n][m] )
{
   int i,j;

   for(i=0; i<n; i++)
      for(j=0; j<m; j++)
         a[i][j] = create1(cos(rand()));
}



void copy_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m] )
{
   int i,j;

   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
         b[i][j] = a[i][j];
}
void multiply_dcmatrix ( int n, int m, int l, dcmplx a[n][m], dcmplx b[m][l], dcmplx c[n][l] )
{
    int i, j, k;
    for ( i=0; i<n; i++ )
	{
      for ( j=0; j<l; j++ )
	  {
           c[i][j] = zero;
           for ( k=0; k<m; k++ )
		  {
		    c[i][j]=add_dcmplx( c[i][j], mul_dcmplx(a[i][k], b[k][j]));
		  }
	  }
	}
}

void add_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m], dcmplx c[n][m] )
{
    int i, j;
    for( i=0; i<n; i++)
      for( j=0; j<m; j++)
        {
          c[i][j]=add_dcmplx(a[i][j], b[i][j]);
        }
}

void sub_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m], dcmplx c[n][m] )
{
    int i, j;
    for( i=0; i<n; i++)
      for( j=0; j<m; j++)
        {
          c[i][j]=sub_dcmplx(a[i][j], b[i][j]);
        }
}

   
void free_dcmatrix ( int n, dcmplx ** a)
{
  int i, j;
  for(i=0; i<n; j++)
       free(a[i]);
  free(a);
}

void swap_dc ( dcmplx *x, dcmplx *y )
{
   dcmplx tmp = *y;
   *y = *x;
   *x = tmp;
}

int lufac ( int n, dcmplx a[n][n], int ipvt[n] )
{
   int info,k,kp1,p,nm1,i,j;
   double smax,dfacc;
   dcmplx dcacc,zero;

   info = -1;
   nm1 = n-1;
   if (nm1 >= 1)
   {
      for (k=0; k<nm1; k++)
      {
         kp1 = k+1;
         p = k;                        /* find pivot index p */
         smax = modulus(a[k][k]);
         for (i=kp1; i<n; i++)
         {
            dfacc = modulus(a[i][k]);
            if (dfacc > smax)
            {
               p = i;
               smax = dfacc;
            }
         }
         ipvt[k] = p;
         if (smax == 0.0)
            info = k;                  /* column is already triangulated */
         else
         {
            if (p != k) swap_dc(&a[p][k],&a[k][k]);
            dcacc = create1(-1.0);     /* compute multipliers */
            dcacc = div_dcmplx(dcacc,a[k][k]);
            for (i=kp1; i<n; i++)
               a[i][k] = mul_dcmplx(a[i][k],dcacc);
            for (j=kp1; j<n; j++)      /* row elimination */
            {
               if (p != k) swap_dc(&a[p][j],&a[k][j]);
               for (i=kp1; i<n; i++)
                  a[i][j] = add_dcmplx(a[i][j],mul_dcmplx(a[k][j],a[i][k]));
            }
	 }
      }
   }
   ipvt[n-1] = n-1;
   if (modulus(a[n-1][n-1]) == 0.0) info = n-1;

   return info;
}

void lusolve ( int n, dcmplx a[n][n], int ipvt[n], dcmplx b[n] )
{
   int k,p,nm1,kb,i,j;
   dcmplx tmp,acc;

   nm1 = n-1;
   if (nm1 >= 1)                     /* solve l*y = b */ 
   {
      for (k=0; k<nm1; k++)
      {
	 p = ipvt[k];
         tmp = b[p];
         if (p != k)
	 {
            b[p] = b[k];
            b[k] = tmp;
         }
         for (i=k+1; i<n; i++)
	    b[i] = add_dcmplx(b[i],mul_dcmplx(tmp,a[i][k]));
      }
   }
   for (k=0; k<n; k++)               /* solve u*x = y */ 
   {
      kb = n - 1 - k;
      b[kb] = div_dcmplx(b[kb],a[kb][kb]);
      tmp = min_dcmplx(b[kb]);
      for (j=0; j<kb; j++)
         b[j] = add_dcmplx(b[j],mul_dcmplx(tmp,a[j][kb]));
   }
}

dcmplx determinant( int n, dcmplx (*a)[n] )
{
  dcmplx lu[n][n], result;
  int ipvt[n], info, i, sign;

  sign = 1;
  copy_dcmatrix(n, n, a, lu);
  info = lufac(n, lu , ipvt);

  if (info != -1)
    return zero;     /* the matrix is singular */ 
  else
  {
     for(i=0; i<n; i++)
     {
       if(i!=ipvt[i])
         sign = sign * (-1);
     }
    result = one;
    for( i=0; i<n; i++)
      result = mul_dcmplx(result, lu[i][i]);
    result = mul_dcmplx(create1(sign), result);
    return result;
  }    

}

void I_dcmatrix(int n, dcmplx a[n][n])
{
  int i, j;
  
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      {
        if(i==j)
          a[i][j]=one;
        else
          a[i][j]=zero;
      }
}
