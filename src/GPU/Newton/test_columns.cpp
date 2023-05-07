/* Tests making the column representation of a polynomial system.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <cstdlib>
#include "cyclic_columns.h"

using namespace std;

int main ( void )
{
   cout << "testing column representation of cyclic n-roots ..." << endl;

   cout << "give the dimension : ";
   int dim; cin >> dim;

   int **nvr = new int*[dim];
   for(int i=0; i<dim; i++) nvr[i] = new int[dim];

   make_cyclic_variables(dim,dim,nvr);

   int ***idx = new int**[dim]; // we have dim columns and
   for(int i=0; i<dim; i++)     // dim monomials in each column
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++) idx[i][j] = new int[nvr[i][j]];
   }
   make_cyclic_columns(dim,dim,nvr,idx);

   write_cyclic_columns(dim,dim,nvr,idx);

   return 0;
}
