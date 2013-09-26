/* file dc_interpolation.c provides an implementation of the Divided Difference method */

#include "dcmplx.h"
#include "dc_interpolation.h"


void divided_difference ( int n, dcmplx x[n], dcmplx f[n], int m[n] )
{
   int i, j, k, t;
   dcmplx t1, t2, ff[n], tf[n], tmp;
   double tol=1.0e-8;
 
   
   for(i=0; i<n; i++)
     tf[i]=f[i];
 
   for(i=0; i<n; i++)
   {  
     if(m[i]>1)
     {
       t=m[i];
       tmp=f[i];
       while(t>1)
       {
         t--;
         f[++i]=tmp;
       }
     } 
  }

   for(j=0; j<n-1; j++)
   {
     for(i=n-1; i>j; i--)
     {  
       if(!equal_dcmplx(x[i], x[i-j-1], tol))
       { 
         t1 = sub_dcmplx (f[i], f[i-1]);
         t2 = sub_dcmplx (x[i], x[i-j-1]);
         f[i] = div_dcmplx (t1, t2);
       }
       else
       {
         t=m[i];
         tmp=tf[i-m[i]+2+j];

         while(t-j>1)
	 {
           t--;
           f[i--]=tmp; 
	 }
         i++;
       }
     }

     /* 
     printf("the %dth iteration:\n", j);
     for(k=0; k<n; k++)
       writeln_dcmplx(f[k]);
     */    
   }
   
   for(k=0; k<n; k++)
     ff[k]=f[k];

   f[0]=ff[n-1];
   for(i=1; i<=n-1; i++)
   { 
     f[i]=f[i-1];
     for(j=i-1; j>0; j--)
     {
       f[j] = min_dcmplx (mul_dcmplx (f[j], x[n-i-1]));
       f[j] = add_dcmplx (f[j], f[j-1]);
     }
    f[0] = min_dcmplx (mul_dcmplx (f[0], x[n-i-1]));
    f[0] = add_dcmplx (f[0], ff[n-i-1]);
   }
   
}



