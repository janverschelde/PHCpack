/* file dc_roots.c provides an implementation of the method of Weierstrass */

#include <stdlib.h>
#include "dc_roots.h"

void roots ( int n, dcmplx p[n], double eps, int maxit, dcmplx z[n-1] )
{
   dcmplx pp[n],dz,dp;
   double mpz,mdz,maxpz,maxdz;
   int i,k;

   for(i=0; i<n-1; i++)
   {
      z[i] = random_dcmplx1();
      pp[i] = div_dcmplx(p[i],p[n-1]);
   }
   pp[n-1] = create1(1.0);

   for(k=0; k<maxit; k++)
   {
      maxpz = 0.0;
      maxdz = 0.0;
      for(i=0; i<n-1; i++)
      {
         dz = horner(n,pp,z[i]);
         /* printf("residual at z[%d] : ", i);
            writeln_dcmplx(dz); */
         mpz = modulus(dz);
         if (mpz > maxpz) maxpz = mpz;
         dp = difference_product(n-1,i,z);
         dz = div_dcmplx(dz,dp);
         mdz = modulus(dz);
         if (mdz > maxdz) maxdz = mdz;
         /* printf("  update to z[%d] : ", i);
            writeln_dcmplx(dz); */
         z[i] = sub_dcmplx(z[i],dz);
      }
      if ((maxpz <= eps) || (maxdz <= eps)) break;
   }
}

void multiplicities ( int n, double tol, dcmplx r[n], int m[n] )
{
   int i,j;
   dcmplx delta;

   for(i=0; i<n; i++) m[i] = 1;

   for(i=0; i<n; i++)
      for(j=i+1; j<n; j++)
      {
         delta = sub_dcmplx(r[i],r[j]);
         if (dcabs(delta) < tol) 
         {  
             m[i]++; 
             m[j]++;
         }
      }
}

void multiple_roots ( int n, dcmplx p[n], double eps, int maxit,
                      dcmplx r[n-1], double tol, int m[n-1] )
{
   int i,nit;
   dcmplx *dp[n-1];

   roots(n,p,eps,maxit,r);
   multiplicities(n-1,tol,r,m);
   derivatives(n,n-1,p,dp);
   /* write_derivatives(n,n-1,dp); */

   for(i=0; i<n-1; i++)
   {
      if (m[i] == 1)
         nit = Newton(n,p,dp[0],&r[i],eps,8);
      else
         nit = Newton(n-m[i]+1,dp[m[i]-2],dp[m[i]-1],&r[i],eps,8);

      /* printf("refined root %d : ", i);
         writeln_dcmplx(r[i]); */
   }
}
   
dcmplx horner ( int n, dcmplx *p, dcmplx x )
{
   dcmplx res = p[n-1];
   int i;

   for(i=n-2; i>= 0; i--)
   {
      res = mul_dcmplx(res,x);
      res = add_dcmplx(res,p[i]);
   }

   return res;
}

int Newton ( int n, dcmplx p[n], dcmplx dp[n-1], dcmplx *z,
             double eps, int maxit )
{
   int i;
   dcmplx delta,y,dy;

   for(i=0; i<maxit; i++)
   {
      y = horner(n,p,*z);
      dy = horner(n-1,dp,*z);
      delta = div_dcmplx(y,dy);
      *z = sub_dcmplx(*z,delta);
      if(dcabs(delta) <= eps) return i;
   }

   return -1;
}


dcmplx difference_product ( int n, int i, dcmplx z[n] )
{
   dcmplx res = create1(1.0);
   int j;

   for(j=0; j<n; j++)
      if (j != i)
         res = mul_dcmplx(res,sub_dcmplx(z[i],z[j]));

   return res;
}

void derivative ( int n, dcmplx p[n], dcmplx dp[n-1] )
{
   int i;
   double df_power;
   dcmplx dc_power;

   for(i=1; i<n; i++)
   {
      df_power = (double) i;
      dc_power = create1(df_power);
      dp[i-1] = mul_dcmplx(p[i],dc_power);
   }
}

void derivatives ( int n, int m, dcmplx p[n], dcmplx *dp[m] )
{
   int i,j,work_n = n;
   dcmplx work_p[n],work_dp[n-1];

   for(i=0; i<n; i++)
      work_p[i] = p[i];

   for(i=0; i<m; i++)
   {
      dp[i] = (dcmplx*) calloc(--work_n,sizeof(dcmplx));
      derivative(work_n+1,work_p,work_dp);
      for(j=0; j<work_n; j++)
      {
         dp[i][j] = work_dp[j];
         work_p[j] = work_dp[j];
      }
   }
}

void write_derivatives ( int n, int m, dcmplx *dp[m] )
{
   int i,j,nn;

   for(i=0,nn=n-1; i<m; i++,nn--)
   {
      printf("derivative %d :\n", i+1);
      for(j=0; j<nn; j++)
         writeln_dcmplx(dp[i][j]);
   }
}
