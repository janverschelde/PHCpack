#include <stdlib.h>
#include "realization.h"
#include "dcmplx.h"
#include "poly_smith.h"
#include "poly_hermite.h"
#include "poly_gcd.h"
#include "poly_matrix.h"
#include "dc_matrix.h"
#include "dc_inverse.h"
#include "append_polymatrix.h"
#include "dc_roots.h" /* for test only */

#define tol 10e-8

void limit(int n, int m, POLY dd, POLY MM[n][m], dcmplx D[n][m] )
{
  int i, j, d_deg, n_deg;
  
  for ( i=0; i<n; i++ )
    for ( j=0; j<m; j++ )
    {
      d_deg = dd.d;
      n_deg = MM[i][j].d;
      if( d_deg < n_deg )
	{
          printf("realizations do not exist\n");
          return;
        }
      if( d_deg == n_deg )
        {
          D[i][j] = div_dcmplx ( MM[i][j].p[n_deg], dd.p[d_deg] );
        }
      if( d_deg > n_deg )
        {
          D[i][j] = create1(0.0);
        }
     }
}


int get_deg(int n, POLY M1[n][n], int a[n], int *sum)
{
  int i, j, flag;

  *sum=0;
  flag=1;
  for( j=0; j<n; j++)
    {
     a[j] = 0;
     for( i=0; i<n; i++)
       {
         if(M1[i][j].d>a[j])
             a[j] = M1[i][j].d;
       }
     if(a[j]==0)
       flag=0;
     *sum = *sum + a[j];
    }
  return flag;
}


void factorize_M1(int p, POLY M1[p][p], int a[p], int n, dcmplx DH[p][p],
         POLY Dia[p][p], dcmplx DL[p][n], POLY S[n][p] )
{
  int i, j, k, col_deg, deg, u, low_bound, up_bound;
    

  for( j=0; j<p; j++ )
     for( i=0; i<p; i++ )
       {
	 /* calculate DH */
         col_deg = a[j];
         deg = M1[i][j].d;
         if( col_deg>deg ) DH[i][j]=create1(0.0);
         if( col_deg==deg) DH[i][j]=M1[i][j].p[deg];
        
         /* get Dia matrix */
         if( i==j )
	   {
             Dia[i][j].d=col_deg;
             Dia[i][j].p=(dcmplx*) calloc(col_deg+1, sizeof(dcmplx));
             for( k=0; k<col_deg; k++ )
                 Dia[i][j].p[k] = zero;   
             Dia[i][j].p[col_deg] = one;
           }
         if( i!=j )
           {
            Dia[i][j].p = (dcmplx*) calloc(1, sizeof(dcmplx));
            Dia[i][j].d = 0;
            Dia[i][j].p[0]=zero;
           }
       }
  /* printf("the DH matrix is:\n"); print_dcmatrix(p, p, DH); */

  /* get the S matrix */
   k=0; deg=0; u=0;
   low_bound;
   up_bound;
   for( j=0; j<p; j++ )
     {
       low_bound = k;
       up_bound = k+a[j];
       for(i=0; i<n; i++)
	 {   
           if(i<up_bound && i>=low_bound)   
	     {
               deg = i - low_bound;
               S[i][j].d = deg;
               S[i][j].p = (dcmplx*) calloc(deg+1, sizeof(dcmplx));
               for(u=0; u<deg; u++)
		  S[i][j].p[u] = zero;
               S[i][j].p[deg] = one; 
               k++;
	     }
           else
	     {
               S[i][j].d = 0;
               S[i][j].p = (dcmplx*) calloc(1, sizeof(dcmplx));
               S[i][j].p[0]=zero;
	     }
	 }
     }
 
   /* printf("the S matrix is:\n"); print1(n, p, S); */

   /* calculate DL matrix */

    for( i=0; i<p; i++ )
      {
        j=0;   k=0;
        low_bound=0;
        while((k<p) && (a[k]==0))
        {
          k++;
	}
	while( j<n && k<p)
	  {
            deg = j-low_bound;   
            if((deg<a[k])&&(deg>=0))
	      {
               
                DL[i][j]=M1[i][k].p[deg]; 
                j++;
                
              }
            else
             {
               DL[i][j] = zero;
               j++;
             }  
            if(((j-low_bound) == a[k] )&& (a[k]!=0))
              {   
                  low_bound =low_bound+a[k];
                  k++;
	      }
	    if((k<p)&&(a[k]==0)) k++;   
	  }
      }       
    /* printf("the DL matrix is:\n"); print_dcmatrix(p, n, DL); */
}
  

POLY Get_Hs (int p, int m,  POLY M1[p][p], POLY M2[m][p], POLY Hs[m][p])
{
  POLY ds;
  dcmplx ds_value, s, H[m][p];
  int i, j;   

  ds = Inverse_Poly (p, M1);

  /* calculate M2 * (M1's inverse matrix) = ds * Hs */
  Multiply(m, p, p, M2, M1, Hs);

  /*********************************************************/
  /***** The following part is for test only ***************/
  /*********************************************************/ 

  /* evaluate the value of the H function at the point s */
  /* s = create2(-0.23423423423, 0); 
     ds_value = horner(ds.d+1, ds.p, s); */ 

  /* initialize the matrix H */
  /* zero_dcmatrix(m, p, H);
 
  /* evaluate_matrix(m, p, Hs, H, s);
     for(i=0; i<m;i++)
       for(j=0; j<p; j++)
         H[i][j] = div_dcmplx(H[i][j],ds_value); */

 /* 
  printf("\nThe H matrix value before realization at the point s is:\n");
  printf("s = ");
  writeln_dcmplx(s);
  print_dcmatrix(m, p, H);   
 */ 
  /*
  printf("H[0][0] =  ");
  write_dcmplx(H[0][0]);
  if(dabs(H[0][0].im-0.0) > 10e-5){ printf("  Complex!!!!!!!!!!\n\n");  getchar(); }
  else   printf("  Real\n\n");
  */
  
  return ds;
}


void Get_B_c( int sum, int p, int d[p], dcmplx B_c[sum][p] )
{
  int i, j, counter;
  counter = 0; 
  for(j=0; j<p; j++)
    {
      counter = counter + d[j];
      for(i=0; i<sum; i++)
       {       
         if((i==(counter-1))&& (d[j]!=0)) B_c[i][j]=one;
         else  B_c[i][j]=zero;
       }
    }
}

void Get_Bc1(int sum, int new_p, int p, int b[new_p],dcmplx DH[new_p][new_p], dcmplx Bc[sum][p])
{
  int i, j;
  dcmplx t_Bc[sum][new_p];

  Get_Bc(sum, new_p, b, DH, t_Bc);

  for(j=0; j<new_p; j++)
    for(i=0; i<sum; i++)
      Bc[i][j]=t_Bc[i][j];  
   
}


void Get_Bc ( int sum, int p, int d[p], dcmplx DH[p][p], dcmplx Bc[sum][p])
{
  int i;
  dcmplx Bm[p][p], B_c[sum][p];

  copy_dcmatrix(p, p, DH, Bm); 
  dcinverse(p, Bm);
  Get_B_c(sum, p, d, B_c);
  /* printf("B_c is:\n"); print_dcmatrix(sum, p, B_c); */
  /* printf("Bm is:\n"); print_dcmatrix(p, p, Bm); */   
  multiply_dcmatrix(sum, p, p, B_c, Bm, Bc);
  /* printf("The matrix Bc is:\n"); print_dcmatrix(sum, p, Bc); */

}

void Get_A_c( int sum, int p, int d[p], dcmplx A_c[sum][sum])
{
  int i, j, k, counter[p];
  counter[0]=d[0];
  for( i=1; i<p; i++)
    counter[i]=counter[i-1]+d[i];
  k=0;
  for( i=0; i<sum; i++)
    {
     if(d[k]==0) k++;
     for( j=0; j<sum; j++)
      { 
        if((i!=(counter[k]-1)) && (j==i+1)&& (d[k]!=1))
           A_c[i][j]=one;
        else
           A_c[i][j]=zero;    
      }
     if(i==(counter[k]-1)) k++;
    }
 
}

void Get_Ac ( int sum, int p, int d[p], dcmplx DH[p][p], dcmplx DL[p][sum], dcmplx Ac[sum][sum] )
{
  dcmplx DH_inverse[p][p], Am[p][sum], neg_Am[p][sum], tmp[sum][sum],
         B_c[sum][p];
  dcmplx A_c[sum][sum]; 
  int i,j;

  copy_dcmatrix(p, p, DH, DH_inverse);
  /* printf("the DH matrix is:\n"); print_dcmatrix(p, p, DH); */
  dcinverse(p, DH_inverse);
  /* printf("the DL matrix is:\n"); print_dcmatrix(p, sum, DL); */
  multiply_dcmatrix(p, p, sum, DH_inverse, DL, neg_Am);
  for(i=0; i<p; i++)
    for(j=0; j<sum; j++)
      Am[i][j]=min_dcmplx(neg_Am[i][j]);

  Get_B_c(sum, p, d, B_c);
  multiply_dcmatrix(sum, p, sum, B_c, Am, tmp);
  Get_A_c(sum, p, d, A_c);
  /* printf("the A_c matrix is:\n"); print_dcmatrix(sum, sum, A_c); */
  /* printf("the Am matrix is:\n"); print_dcmatrix(p, sum, Am); */
  add_dcmatrix(sum, sum, A_c, tmp, Ac);
  /* printf("The matrix Ac is:\n"); print_dcmatrix(sum, sum, Ac); */
}  

void Get_Cc ( int sum, int m, int p, POLY N[m][p], POLY D[p][p], int d[p], dcmplx Dc[m][p], dcmplx Cc[m][sum])
{
   POLY Dc_D[m][p];
   POLY tmp[m][p];
   int i, j, k, low_bound, deg;
 
   dcmatrix_Multiply ( m, p, p, Dc, D, Dc_D );
   sub_polymatrix(m, p, N, Dc_D, tmp);
   /* printf("the Cc*S=N-Dc*D matrix is:\n"); print1(m, p, tmp); */

    for( i=0; i<m; i++ )
      {
        k=0; j=0;
        low_bound=0;
        while((k<p) && (d[k]==0)) k++;
	while( j<sum && k<p)
	  {
            deg = j-low_bound;
            if((deg<=tmp[i][k].d)&& (deg>=0))
	      {               
                Cc[i][j]=tmp[i][k].p[deg]; 
                j++;                
              }
            else
             {
               Cc[i][j] = zero;
               j++;
             }  
            if(( (j-low_bound) == d[k] )&& (d[k]!=0))
              {   
                  low_bound =low_bound+d[k];
                  k++;
	      }
            if((k<p) && (d[k]==0)) k++;
	  }
      }       
    /* printf("The matrix Cc is:\n"); print_dcmatrix(m, sum, Cc); */
    free_matrix(m, p, Dc_D);
    free_matrix(m, p, tmp);
 
} 
void Get_Cc1 ( int sum, int m, int p, POLY N[m][p], int d[p], dcmplx Cc[m][sum])
{

   int i, j, k, low_bound, deg;

    for( i=0; i<m; i++ )
      {
        k=0; j=0;
        low_bound=0;
	while( j<sum && k<p)
	  {
            deg = j-low_bound;
            
            if(deg<=N[i][k].d)
	      {
               
                Cc[i][j]=N[i][k].p[deg]; 
                j++;
                
              }
            else
             {
               Cc[i][j] = zero;
               j++;
             }  
            if( (j-low_bound) == d[k] )
              {   
                  low_bound =low_bound+d[k];
                  k++;
	      }
	  }
      }       

 
} 

void realization(int p, int m, int q, POLY M1[p][p], POLY M2[m][p], dcmplx Ac[q][q], dcmplx Bc[q][p],
                 dcmplx Cc[m][q], dcmplx Dc[m][p])
{
   POLY    inverse_M1[p][p], Dia[p][p], H[m][p], S[q][p], ds;
   dcmplx  DH[p][p], DL[p][q];
   int     a[p], b[p], c[p], sum0, sum, i, j, new_p, flag;
   POLY**  new_Hs;
   POLY N[m][p], D[p][p], Gr[p][p], tmp, tmp1;
   double  max_mod;

   POLY ds_Dc[m][p], t_ds_Dc[p][p], C_B[m][p];
   dcmplx m1, m2;

   dcmplx t_C_B[m][p];
  
   /*** get transfer function from PHC output 1/ds * H = M2 * M1_inverse */
   copy( p, p, M1, inverse_M1);
   ds = Get_Hs(p, m, inverse_M1, M2, H);
   printf("The denominator of the transfer function of the dynamic compensator is:\n"); 
   Print_Poly(ds.d, ds.p);
   printf("The numerator (%d*%d matrix) of the transfer function of the dynamic compensator is:\n", m, p); 
   print1(m, p, H);
   limit(m, p, ds, H, Dc);
              
   get_deg( p, M1, a, &sum );
 
   if(sum != q )  /* also need to check if they have the common gcd for general purpose!!!!! */ 
   {              /* this should not be the case from the PHC output */
     factorize_H1( m, p, ds, H, N, D);
     right_coprime(p, m, p, D, N, Gr);     /* use Hermite form to make N and D matrices right coprime */
     get_deg( p, D, b, &sum );
     factorize_M1( p, D, b, sum, DH, Dia, DL, S );
     Get_Ac(sum, p, b, DH, DL, Ac);
     Get_Bc(sum, p, b, DH, Bc);
     Get_Cc ( sum, m, p, N, D, b, Dc, Cc);
     free_matrix(m, p, N);
     free_matrix(p, p, D);
     return;
   }

  else    /* there may have a zero column degree in M1, we can use modified algorithm (from Linear Systems p. 409) to get the minimal realization */
  {
   factorize_M1( p, M1, a, sum, DH, Dia, DL, S );
   Get_Ac(sum, p, a, DH, DL, Ac);
   Get_Bc(sum, p, a, DH, Bc);
   Get_Cc (sum, m, p, M2, M1, a, Dc, Cc);
  }

  free(ds.p);
  free_matrix(p, p, inverse_M1);
  free_matrix(m, p, H);
  free_matrix(p, p, Dia);
  free_matrix(q, p, S);
}

POLY** eliminate_columns(int m, int p, POLY H_s[m][p], int a[p], int new_p)
{
  int i, j, k;
  POLY** new_Hs;
  
  new_Hs = (POLY**) calloc(m, sizeof(POLY*));  
  for(i = 0; i<m; i++)
    new_Hs[i] = (POLY*) calloc(new_p, sizeof(POLY));  

  k=0;
  for(j=0; j<p; j++)
    {
      if(a[j]!=0)
	{
          for(i=0; i<m; i++)
	    {
              new_Hs[i][k] = H_s[i][j];   
	    }
          k++;
         }
    }
  return new_Hs;

}

void factorize_H1(int m, int p, POLY ds, POLY H[m][p], POLY N[m][p], POLY D[p][p])
{
  int i, j;
  POLY ** Hs;
 
  Hs = (POLY**) calloc(m, sizeof(POLY*));  
  for(i = 0; i<m; i++)
    Hs[i] = (POLY*) calloc(p, sizeof(POLY)); 

  for(i=0; i<m; i++)
    for(j=0; j<p; j++)
      Hs[i][j]=H[i][j];
  factorize_H(m, p, ds, Hs, N, D); 
}

void factorize_H(int m, int p, POLY ds, POLY** H, POLY N[m][p], POLY D[p][p])
{
  int i, j, dgcd, dk, dl, dbd, dad;
  POLY l[p], k[p], tmp;
  dcmplx **gcd_coeff1, **gcd_coeff2;


  for(j=0; j<p; j++)
  {  

    /* find the gcd of a column */
    tmp = assign_poly(H[0][j]);
    k[j]= assign_poly(H[0][j]);
    for(i=1; i<m; i++)
    {      
      free(k[j].p);       
      k[j]= get_gcd(tmp, H[i][j]);
      free(tmp.p);
      tmp = assign_poly(k[j]);     
    }

    if((k[j].d==0) && equal_dcmplx(k[j].p[0], one, tol))  
    {
      for(i=0; i<m; i++)
        N[i][j] =assign_poly(H[i][j]);  
      l[j] = assign_poly(ds);
    }
    else
    {
         gcd_coeff1=ExtPolyGcd( ds.d+1, ds.p, k[j].d+1, k[j].p, &dgcd, &dk, &dl, &dbd, &dad);
         l[j].d = dad;
         l[j].p = assign(dad, gcd_coeff1[4]);
         for(i=0; i<m; i++)
         {
           gcd_coeff2 = ExtPolyGcd( H[i][j].d+1, H[i][j].p, dgcd+1, gcd_coeff1[2], &dgcd, &dk, &dl, &dbd, &dad);
           N[i][j].d = dad; 
           N[i][j].p = assign(dad, gcd_coeff2[4]);
           free_gcd_coeff( gcd_coeff2);
         }

         free_gcd_coeff( gcd_coeff1);
    }
  }

  /*
  for(i=0; i<p; i++)
  { 
    printf("the %dth column of the D matrix is :\n", i);
    printf("l[%d]=", i);
    Print_Poly(l[i].d, l[i].p);
  }
  */     
  for(i=0; i<p; i++)
    for(j=0; j<p; j++)
    {
      if(i==j)
        D[i][j]=assign_poly(l[j]);
      else
	{
          D[i][j].d=0;
          D[i][j].p=(dcmplx*) calloc(1, sizeof(dcmplx));
          D[i][j].p[0]=zero;
        }
     }
  
} 



void add_zero_columns(int q, int p, int a[p], dcmplx Bc[q][p])
{
  int i, j, k;
   
  for(j=0; j<p; j++)
  {
    if( a[j] == 0 )
    { 
      if( j==p-1 )
      {
         for(i=0; i<q; i++)
	   Bc[i][j] = zero;
      }
      else
      {
	for(k=p-1; k>=j+1; k--)
	  for(i=0; i<q; i++)
	    Bc[i][k] = Bc[i][k-1];
        for(i=0; i<q; i++)
          Bc[i][k] = zero;
      }
      
        
    } 
  } 

}

void right_coprime(int n1, int n2, int m, POLY p1[n1][m], POLY p2[n2][m], POLY Gr[m][m])
{
  int  n = n1 + n2, i, j;  
  POLY t_Gr[n][m],inv_Gr[m][m], p[n][n], ds, D[n1][m], N[n2][m], *tmp;
  POLY ds1, pp1[n][n], q1[m][m], t_Gr1[n][m], inv_smith1[m][m], inv_Gr1[m][m];


  poly_v_append ( n1, n2, m, p1, p2, t_Gr);

  /* poly_v_append ( n1, n2, m, p1, p2, t_Gr1); */
  /* Smith(n, m, t_Gr1, pp1, q1); */
 
   Hermite(n, m, t_Gr, p);
   /* printf("The Hermite form is:\n"); print(n, m, t_Gr); */
  
  for(i=0; i<m; i++)
   for(j=0; j<m; j++)
   {
      Gr[i][j]=assign_poly(t_Gr[i][j]);
      inv_Gr[i][j]=assign_poly(t_Gr[i][j]);
      /* inv_smith1[i][j]=assign_poly(t_Gr1[i][j]); */
   }
  /* printf("The Gr is:\n"); print(m, m, Gr); */
  
  ds = Inverse_Poly(m, inv_Gr);
  /* printf("ds is:\n"); Print_Poly(ds.d, ds.p); */
  /* printf("The inverse of Gr from Hermite is:\n");  print(m, m, inv_Gr); */  

  /*
  ds1 = Inverse_Poly(m, inv_smith1);
  Multiply(m, m, m, q1, inv_smith1, inv_Gr1);
  printf("ds1 is:\n"); Print_Poly(ds1.d, ds1.p);
  printf("The inverse of Gr from Smith is:\n");  print(m, m, inv_Gr1);
  */

  copy(n1, m, p1, D);
  copy(n2, m, p2, N);
  free_matrix(n1, m, p1);
  free_matrix(n2, m, p2);
  Multiply(n1, m, m, D, inv_Gr, p1);
  Multiply(n2, m, m, N, inv_Gr, p2);

  /* printf("p1 is:\n"); print(n1, m, p1); */
  /* printf("p2 is:\n"); print(n2, m, p2);  */

  /*check if they are coprime with smith form */
  /*  poly_v_append ( n1, n2, m, p1, p2, t_Gr1); */
  /* Smith(n, m, t_Gr1, pp1, q1); */ 
  /* printf("The Smith form is:\n"); print(n, m, t_Gr1); */
 

  /* printf("\np1 afer coprime:\n"); */
  for(i=0; i<n1; i++)
   for(j=0; j<m; j++)
   {
     tmp=ExtPolyGcd1(ds, p1[i][j]);
     free(p1[i][j].p);
     p1[i][j]=assign_poly(tmp[3]);
     /* printf("p1%d%d=\n",i,j); Print_Poly(tmp[3].d, tmp[3].p); */
     free_gcd_coeff1(tmp);
   }
  for(i=0; i<n2; i++)
   for(j=0; j<m; j++)
   {
     tmp=ExtPolyGcd1(ds, p2[i][j]);
     free(p2[i][j].p);
     p2[i][j]=assign_poly(tmp[3]); 
     free_gcd_coeff1(tmp);     
   }

  free(ds.p);
  free_matrix(n1, m, D);
  free_matrix(n2, m, N);
  free_matrix(m, m, inv_Gr);
  free_matrix(n, m, t_Gr);
  free_matrix(n, n, p);   
}




