// defines code for execution on the host

#include <cmath>
#include "gram_host.h"
#include "gqd_qd_utilT.h"

void print_gram_matrices ( complexH<T1> **g, complex<T> *g_h, int dim )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         complexH<T1> temp;
         comp1_gqd2qd(&g_h[i*dim+j],&temp);
         cout << "GPU g[" << i << "," << j << "] = " << temp;
         cout << "CPU g[" << i << "," << j << "] = " << g[i][j];
      }
}

void print_difference ( complexH<T1> **g, complex<T> *g_h, int dim )
{
   complexH<T1> sum(0.0,0.0);

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         complexH<T1> temp;
         comp1_gqd2qd(&g_h[i*dim+j],&temp);
         temp = temp - g[i][j];
         sum = sum + temp.adj()*temp;
      }
   cout << "2-norm of componentwise difference : " << sqrt(sum.real) << endl;
}

void CPU_inner_product
 ( complexH<T1>** v, int dim, int i, int j, complexH<T1>* ip )
{
   complexH<T1> sum(0.0,0.0);

   for(int k=0; k<dim; k++)
      sum = sum + v[i][k].adj()*v[j][k];

   *ip = sum;
}

void CPU_gram ( complexH<T1>** v, complexH<T1>** g, int dim )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         CPU_inner_product(v,dim,i,j,&g[i][j]);
}
