// Simple test program to evaluate and differentiate 
// the polynomial system chandra for visual inspection
// and verification with a computer algebra system such as Maple.

#include "chandra.h"

using namespace std;

int main ( void )
{
   cout << "give the dimension : ";
   int dim; cin >> dim;

   complexH<T1>* x = new complexH<T1>[dim];
   complexH<T1>* y = new complexH<T1>[dim];
   complexH<T1> c_a;
   c_a.init(33.0,0.0);
   complexH<T1> c_b;
   c_b.init(64.0,0.0);
   complexH<T1> c = c_a/c_b;  // avoid representation errors, unlike 0.51234

   for(int i=0; i<dim; i++) x[i].init(1.0,0.0);
   cout << "the values for x :" << endl;
   for(int i=0; i<dim; i++)
      cout << "->  x[" << i << "] = " << x[i];
   chandra_evaluate(c,dim,x,y);
   cout << "the evaluated x :" << endl;
   for(int i=0; i<dim; i++)
      cout << "->  y[" << i << "] = " << y[i];
   cout << "diagonal elements on the Jacobian :" << endl;
   for(int i=0; i<dim; i++)
   {
      chandra_differentiate(c,dim,i,i,x,&y[0]);
      cout << "-> jf[" << i << "," << i << "] = " << y[0];
   }
   cout << "off diagonal elements on the Jacobian : " << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         if(i != j)
         {
            chandra_differentiate(c,dim,i,j,x,&y[0]);
            cout << "-> jf[" << i << "," << j << "] = " << y[0];
         }
   complexH<T1>** v = new complexH<T1>*[dim+1];
   for(int i=0; i<dim+1; i++) v[i] = new complexH<T1>[dim];
   chandra_evaluate_and_differentiate(c,dim,x,v);
   cout << "The Jacobian matrix : " << endl;
   for(int i=0; i<dim; i++) 
      for(int j=0; j<dim; j++) 
         cout << "-> jf[" << i << "," << j << "] = " << v[j][i];
   cout << "the evaluated system with flipped sign : " << endl;
   for(int i=0; i<dim; i++) 
      cout << "-> f[" << i << "] = " << v[dim][i];

   return 0;
}
