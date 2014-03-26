// Definitions of the prototypes declared in chandra.h.

#include "chandra.h"

void complex_chandra_evaluate
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

void real_chandra_evaluate ( T1 c, int dim, int i, T1* x, T1* y )
{
   T1 onedim,two,twodim,sum;

   two = 2.0;
   onedim = dim;
   twodim = two*onedim;

   *y = twodim*x[i];
   sum = 1.0;
   for(int j=0; j<dim-1; j++)
   {
      T1 cff;
      double d_i = (double) (i+1);
      double d_ipj = (double) (i+1+j+1);
      T1 c_a = d_i;
      T1 c_b = d_ipj;
      cff = c_a/c_b;
      sum = sum + cff*x[j];
   } 
   *y = *y - c*x[i]*sum;
   *y = *y - twodim;
}

void complex_chandra_evaluate
 ( complexH<T1> c, int dim, complexH<T1>* x, complexH<T1>* y )
{
   for(int i=0; i<dim; i++) complex_chandra_evaluate(c,dim,i,x,&y[i]);
}

void real_chandra_evaluate ( T1 c, int dim, T1* x, T1* y )
{
   for(int i=0; i<dim; i++) real_chandra_evaluate(c,dim,i,x,&y[i]);
}

void complex_chandra_differentiate
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

void real_chandra_differentiate ( T1 c, int dim, int i, int j, T1* x, T1* y )
{
   T1 onedim,two,twodim,sum,cff;

   two = 2.0;
   onedim = dim;
   twodim = two*onedim;

   if(i == j)
   {
      sum = 1.0;
      for(int k=0; k<dim-1; k++)
      {
         double d_i = (double) (i+1);
         double d_ipk = (double) (i+1+k+1);
         T1 c_a = d_i;
         T1 c_b = d_ipk;
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
         (*y) = 0.0;
      else
      {
         double d_i = (double) (-i-1);
         double d_ipj = (double) (i+1+j+1);
         T1 c_a = d_i;
         T1 c_b = d_ipj;
         cff = c_a/c_b;
         *y = cff*c*x[i];
      }
   }
}

void complex_chandra_evaluate_and_differentiate
 ( complexH<T1> c, int dim, complexH<T1>* x, complexH<T1>** v )
{
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
         complex_chandra_differentiate(c,dim,i,j,x,&v[j][i]);

   complexH<T1> zero;
   zero.init(0.0,0.0);
   for(int i=0; i<dim; i++)
   {
      complex_chandra_evaluate(c,dim,i,x,&v[dim][i]);
      v[dim][i] = zero - v[dim][i];
   }
}

void real_chandra_evaluate_and_differentiate ( T1 c, int dim, T1* x, T1** v )
{
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
         real_chandra_differentiate(c,dim,i,j,x,&v[j][i]);

   T1 zero = 0.0;
   for(int i=0; i<dim; i++)
   {
      real_chandra_evaluate(c,dim,i,x,&v[dim][i]);
      v[dim][i] = zero - v[dim][i];
   }
}
