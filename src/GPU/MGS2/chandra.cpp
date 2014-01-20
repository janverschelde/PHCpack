// Definitions of the prototypes declared in chandra.h.

#include "chandra.h"

void chandra_evaluate
 ( complexH<T1> c, int dim, int i, complexH<T1>* x, complexH<T1>* y )
{
   complexH<T1> onedim,two,twodim,sum;

   two.init(2.0,0.0);
   onedim.init(dim,0.0);
   twodim = two*onedim;

   *y = twodim*x[i];
   sum.init(1.0,0.0);
   for(int j=0; j<dim-1; j++)
   {
      complexH<T1> cff;
      double d_i = (double) (i+1);
      double d_ipj = (double) (i+1+j+1);
      T1 a = d_i;
      T1 b = d_ipj;
      complexH<T1> c_a;
      complexH<T1> c_b;
      c_a.init(d_i,0.0);
      c_b.init(d_ipj,0.0);
      cff = c_a/c_b;
      sum = sum + cff*x[j];
   } 
   *y = *y - c*x[i]*sum;
   *y = *y - twodim;
}

void chandra_evaluate
 ( complexH<T1> c, int dim, complexH<T1>* x, complexH<T1>* y )
{
   for(int i=0; i<dim; i++) chandra_evaluate(c,dim,i,x,&y[i]);
}

void chandra_differentiate
 ( complexH<T1> c, int dim, int i, int j, complexH<T1>* x, complexH<T1>* y )
{
   complexH<T1> onedim,two,twodim,sum,cff;

   two.init(2.0,0.0);
   onedim.init(dim,0.0);
   twodim = two*onedim;

   if(i == j)
   {
      sum.init(1.0,0.0);
      for(int k=0; k<dim-1; k++)
      {
         double d_i = (double) (i+1);
         double d_ipk = (double) (i+1+k+1);
         T1 a = d_i;
         T1 b = d_ipk;
         complexH<T1> c_a;
         complexH<T1> c_b;
         c_a.init(d_i,0.0);
         c_b.init(d_ipk,0.0);
         cff = c_a/c_b;
         sum = sum + cff*x[k];
      }
      *y = twodim - c*sum;
      if(i < dim-1)
         *y = *y - (c/two)*x[i];
   }
   else
   {
      if(j == dim-1)
         (*y).init(0.0,0.0);
      else
      {
         double d_i = (double) (-i-1);
         double d_ipj = (double) (i+1+j+1);
         T1 a = d_i;
         T1 b = d_ipj;
         complexH<T1> c_a;
         complexH<T1> c_b;
         c_a.init(d_i,0.0);
         c_b.init(d_ipj,0.0);
         cff = c_a/c_b;
         *y = cff*c*x[i];
      }
   }
}

void chandra_evaluate_and_differentiate
 ( complexH<T1> c, int dim, complexH<T1>*x, complexH<T1>** v )
{
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
         chandra_differentiate(c,dim,i,j,x,&v[j][i]);

   complexH<T1> zero;
   zero.init(0.0,0.0);
   for(int i=0; i<dim; i++)
   {
      chandra_evaluate(c,dim,i,x,&v[dim][i]);
      v[dim][i] = zero - v[dim][i];
   }
}
