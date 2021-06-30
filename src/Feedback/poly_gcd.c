/* poly_gcd.c provides an implementation of the extended gcd method */

#include <stdlib.h>
#include "dc_roots.h"
#include "poly_gcd.h"
#include "poly_dcmplx.h"
#include "dc_interpolation.h"
#include "dcmplx.h"

POLY* ExtPolyGcd1 ( POLY a, POLY b )
{
   dcmplx **gcd_coeff;
   POLY* result=(POLY*) calloc(5, sizeof(POLY));

   int i, dgcd, dk, dl, dbd, dad;

   gcd_coeff=ExtPolyGcd(a.d+1, a.p, b.d+1, b.p, &dgcd, &dk, &dl, &dbd, &dad);

   for(i=0; i<5; i++) result[i].p=gcd_coeff[i];

   result[0].d=dk;
   result[1].d=dl;
   result[2].d=dgcd;
   result[3].d=dbd;
   result[4].d=dad;
  
   return result;
}

dcmplx** ExtPolyGcd
 ( int n, dcmplx a[n], int m, dcmplx b[m], int *dgcd, int *dk, int *dl,
   int *dbd, int *dad )
{
   dcmplx **result=(dcmplx**) calloc(5, sizeof(dcmplx*));  
   dcmplx ra[n-1], rb[m-1];
   dcmplx *k, *l,*rbd, *rad, *ad, *bd, *g_cd; 
   dcmplx coef_aa, coef_bb;
   int i, j;
   double tol=1.0e-8;

   if((n==1)||(m==1))
   { 
    /* suppose the gcd of any polynomial and 0 is the polynomial itself */
       if(equal_dcmplx(b[0], create1(0), tol)&&(m==1))
       {
       k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
       k[0]=create1(1.0);
       l=(dcmplx*) calloc( 1, sizeof(dcmplx) ); 
       l[0]=create1(1.0);

         *dgcd=n-1;
         g_cd=(dcmplx*) calloc( n, sizeof(dcmplx) );
         for(i=0; i<n; i++)
           g_cd[i]=a[i];

         *dbd=0;
         bd=(dcmplx*) calloc(1, sizeof(dcmplx));
         bd[0]=create1(0.0);
         *dad=0;
         ad=(dcmplx*) calloc(1, sizeof(dcmplx));
         ad[0]=create1(1.0); 
       }
       else if(equal_dcmplx(a[0], create1(0.0), tol)&&(n==1))
       {
       k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
       k[0]=create1(1.0);
       l=(dcmplx*) calloc( 1, sizeof(dcmplx) ); 
       l[0]=create1(1.0);

         *dgcd=m-1;
         g_cd=(dcmplx*) calloc( m, sizeof(dcmplx) );
         for(i=0; i<m; i++)
           g_cd[i]=b[i];
 
         *dbd=0;
         bd=(dcmplx*) calloc(1, sizeof(dcmplx));
         bd[0]=create1(1.0);
         *dad=0;
         ad=(dcmplx*) calloc(1, sizeof(dcmplx));
         ad[0]=create1(0.0); 

      }
    
     else
     {
      if(n==1)
      {    
        k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        k[0]=div_dcmplx(create1(1.0), a[0]);
        l=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        l[0]=create1(0);
      }
      else
      {  
        l=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        l[0]=div_dcmplx(create1(1.0), b[0]);
        k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        k[0]=create1(0); 
      }
     *dgcd=0;
     g_cd=(dcmplx*) calloc( 1, sizeof(dcmplx) );
     g_cd[0]=create1(1.0);
     *dbd=m-1;
     bd=assign(*dbd, b);
     *dad=n-1;
     ad=assign(*dad, a);
     }
     *dk=0;
     *dl=0;

     result[0]=k;
     result[1]=l;
     result[2]=g_cd;
     result[3]=bd;
     result[4]=ad;

     return result;
   }

   rootsGCD(n,a,m,b,dgcd,ra,rb); 
 
    *dk = m-2-*dgcd;
    *dl = n-2-*dgcd;

    if(*dk<0 || *dl<0)
    {
      if(*dk<0)
      { 
        *dk=0; *dl=0;
        k=(dcmplx*) calloc (1, sizeof(dcmplx));
        k[0]=create1(0);
        l=(dcmplx*) calloc (1, sizeof(dcmplx));
        l[0]=div_dcmplx(create1(1.0), b[m-1]);
      }
      else
      {
        *dl=0; *dk=0;
        l=(dcmplx*) calloc (1, sizeof(dcmplx));
        l[0]=create1(0);
        k=(dcmplx*) calloc (1, sizeof(dcmplx));
        k[0]=div_dcmplx(create1(1.0), a[n-1]);
      }
    }
    else  
    {
       k=(dcmplx*) calloc( m-1-(*dgcd), sizeof(dcmplx));
       l=(dcmplx*) calloc( n-1-(*dgcd), sizeof(dcmplx));
       extended_gcd ( n-1, ra, m-1, rb, *dgcd, k, l);

      for(i=0; i<=*dk; i++)
        k[i]=div_dcmplx(k[i], a[n-1]);

      for(i=0; i<=*dl; i++)
        l[i]=div_dcmplx(l[i], b[m-1]);  
    }  

      
    *dbd=m-1-*dgcd;
    *dad=n-1-*dgcd;
    rbd=(dcmplx*) calloc(*dbd+1, sizeof(dcmplx));
    rad=(dcmplx*) calloc(*dad+1, sizeof(dcmplx));    
    for(i=0; i<*dbd; i++)
      rbd[i]=rb[i+*dgcd]; 
    for(j=0; j<*dad; j++)
      rad[j]=ra[j+*dgcd];
   
    bd=get_poly(*dbd, rbd);
    ad=get_poly(*dad, rad);
    for(i=0; i<=*dbd; i++)
       bd[i]=mul_dcmplx(b[m-1], bd[i]);

    for(j=0; j<=*dad; j++)
       ad[j]=mul_dcmplx(a[n-1], ad[j]);

    g_cd=get_poly(*dgcd, ra);
   
    if(equal_poly(*dgcd, g_cd, n-1, a)||equal_poly(*dgcd, g_cd, m-1, b))
    {
      free(k); free(l); 
      *dk=0; *dl=0;
      if(equal_poly(*dgcd, g_cd, n-1, a))
      {
         k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
         k[0]=create1(1.0);
         l=(dcmplx*) calloc( 1, sizeof(dcmplx) );
         l[0]=create1(0);
      }
      else
      {
        l=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        l[0]=create1(1.0);
        k=(dcmplx*) calloc( 1, sizeof(dcmplx) );
        k[0]=create1(0);
      }
    }
   
     result[0]=k;
     result[1]=l;
     result[2]=g_cd;
     result[3]=bd;
     result[4]=ad;

   free(rbd); free(rad); 
   return result; 

}

void rootsGCD
 ( int n, dcmplx a[n], int m, dcmplx b[m], int *l,
   dcmplx ra[n-1], dcmplx rb[m-1] )
{
   double eps = 1.0e-12;
   double tol = 1.0e-8;
   double mul_tol = 1.0e-2;

   int i, j, mult1[n-1], mult2[m-1] ;
   int c=0;
  
   multiple_roots(n, a, eps, 10*n, ra, mul_tol, mult1);
   multiple_roots(m, b, eps, 10*m, rb, mul_tol, mult2);

   for(i=c; i<n-1; i++)
      for(j=c; j<m-1; j++)
      {
         if(equal_dcmplx(ra[i], rb[j], tol))
         { 
            if( c!=i ) swap(&ra[c], &ra[i]);
            if( c!=j ) swap(&rb[c], &rb[j]);
            c++;
            break;           
         }
      }

   *l=c;
}

POLY get_gcd ( POLY a, POLY b )
{
   POLY c;
   int dgcd;
   dcmplx ra[a.d], rb[b.d];
  
  /* assume gcd of any polynomial and zero is the polynomial itself */
   if((a.d==0) || (b.d==0))
   {
      if((a.d==0)&&equal_dcmplx(a.p[0], zero, 10e-8))
      {
         c = assign_poly(b);
      }
      else if((b.d==0)&&equal_dcmplx(b.p[0], zero, 10e-8))
      {
         c = assign_poly(a);
      }
      else
      {
         c.d = 0;
         c.p = (dcmplx*) calloc(1, sizeof(dcmplx));
         c.p[0] = create1(1.0);
      }
   }
   else
   {
      rootsGCD(a.d+1,a.p,b.d+1,b.p,&dgcd,ra,rb);
      c.d = dgcd;
      c.p = get_poly(dgcd, ra);
   }
   return c;
}

void extended_gcd
 ( int n, dcmplx ra[n], int m, dcmplx rb[m], int c,
   dcmplx k[m-c], dcmplx l[n-c])
{
   dcmplx kx[m-c], lx[n-c];  
   int i, j, m1[m-c], m2[n-c];
   dcmplx one=create1( 1.0 ); 
   double tol=1.0e-8;

   for(i=0; i<n-c; i++)
   {
      lx[i] = ra[i+c];
   }
   
   for(i=0; i<m-c; i++)
   {
      kx[i] = rb[i+c];
   }

   for(i=0; i<n-c; i++)
   { 
      l[i] = sub_dcmplx( lx[i], kx[0] );
      for(j=1; j<m-c; j++) 
         l[i] = mul_dcmplx( l[i], sub_dcmplx( lx[i], kx[j] ));
      l[i] = div_dcmplx( one, l[i] );
   }

   for(i=0; i<m-c; i++)
   {  
      k[i] = sub_dcmplx( kx[i], lx[0]);
      for(j=1; j<n-c; j++)
         k[i] = mul_dcmplx( k[i], sub_dcmplx( kx[i], lx[j] ));
      k[i] = div_dcmplx( one, k[i] ); 
   }

   if(group_points(n-c, tol, lx, l, m2))
   { 
      hermite_derivative(m-c, kx, n-c, lx, l, m2);
   }  
   if(group_points(m-c, tol, kx, k, m1))
      hermite_derivative(n-c, lx, m-c, kx, k, m1); 

   /*
   printf("the points after group\n");
   for(i=0; i<m; i++)
     writeln_dcmplx(kx[i]);
   */
   divided_difference( n-c, lx, l, m2);
   divided_difference( m-c, kx, k, m1);
} 
     
dcmplx *get_poly ( int n, dcmplx * root )
{ 
   int i, j;
   dcmplx *root1;
   dcmplx *poly=(dcmplx*) calloc(n+1, sizeof(dcmplx));

   root1=assign(n, root);
   if(n==0) 
   {
      poly[0]=create1(1.0);
      return poly;
   }
   poly[0]=min_dcmplx(root1[0]);
   poly[1]=create1(1.0);

   for(i=1; i<n; i++)
   {
      root1[i]=min_dcmplx(root1[i]);
      for(j=i; j>0; j--)
      {
         poly[j]=add_dcmplx(mul_dcmplx(root1[i], poly[j]), poly[j-1]);
      }
      poly[0]=mul_dcmplx(root1[i], poly[0]);
      poly[i+1]=create1(1.0);
   }
   return poly;
}
  
void free_gcd_coeff (dcmplx **gcd_coeff )
{
   int i;
 
   for(i=0; i<5; i++)
   {
      free(gcd_coeff[i]);
   }
   free(gcd_coeff);
}

void free_gcd_coeff1 ( POLY *result )
{
   int i;
 
   for(i=0; i<5; i++)
   {
      free(result[i].p);
   }
   free(result);
}

int group_points ( int n, double tol, dcmplx x[n], dcmplx f[n], int m[n] )
{
   int i, j, tmp, mult;
   dcmplx delta;

   for(i=0; i<n; i++) m[i] = 1;

   mult = 0;
   for(i=0; i<n; i++)
   {
      for(j=i+1; j<n; j++)
      {
         delta = sub_dcmplx(x[i],x[j]);
         if(dcabs(delta) < tol) 
         {  
            m[i]++; 
            m[j]++;
	    if(j!=i+1)
	    {
	       swap(&x[i+1], &x[j]);
               swap(&f[i+1], &f[j]);
               tmp=m[j];
               m[j]=m[i+1];
               m[i+1]=tmp;     
	    }
         }
      }
      if(m[i]>1) mult = 1;
   }
   return mult;
}

dcmplx **rational_derivative
 ( int n, dcmplx *num, int m, dcmplx *denom, int *d_num, int *d_denom )
{  
   dcmplx *t1, *t2;
   int t1_d, t2_d;
   dcmplx dnum[n], ddenom[m];
   dcmplx **derive, *deriv_num, *deriv_denom;
 
   derivative(n+1, num, dnum);  
   derivative(m+1, denom, ddenom);

   t1=mul_poly(n-1, dnum, m, denom, &t1_d);
   t2=mul_poly(n, num, m-1, ddenom, &t2_d);
   deriv_num=min_poly(t1_d, t1, t2_d, t2, d_num);
   deriv_denom=mul_poly(m, denom, m, denom, d_denom);
 
   derive=(dcmplx**) calloc(2, sizeof(dcmplx*));
   derive[0]=deriv_num;
   derive[1]=deriv_denom;
  
   return derive;
}

void hermite_derivative
 ( int root_num, dcmplx root[root_num], int p_num,
   dcmplx x[p_num], dcmplx f[p_num], int m[p_num])
{
   int i,j;
   dcmplx *pp, *num, *denom, n_value, d_value, result;
   dcmplx **deriv;
   int d_num, d_denom, td1, td2, fact; 

   pp=get_poly(root_num, root);
   for(i=0; i<p_num; i++)
   {
      fact=1;
      td1=0;
      td2=root_num;
      if(m[i]>1)
      {       
         num=(dcmplx*) calloc(1, sizeof(dcmplx));
         num[0]=one;
         denom=assign(td2, pp);
         for(j=0; j<m[i]-1; j++)
         {
            fact=fact*(j+1);
            deriv=rational_derivative(td1,num,td2,denom,&d_num,&d_denom);
            n_value=horner(d_num+1, deriv[0], x[i]); 
            d_value=horner(d_denom+1, deriv[1], x[i]);
            result = div_dcmplx(n_value, d_value);
               
            f[++i]=div_double(result,fact);
            free(num);
            free(denom);    
            if(j!=m[i]-2)
            {
               td1=d_num; 
               td2=d_denom;
               num=assign(d_num, deriv[0]);
               denom=assign(d_denom, deriv[1]);
            }
            free(deriv[0]);
            free(deriv[1]);
            free(deriv);
         }
      }
   }
   free(pp); 
}
