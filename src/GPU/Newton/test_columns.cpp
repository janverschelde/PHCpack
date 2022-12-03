/* Tests making the column representation of a polynomial system.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <cstdlib>

using namespace std;

void make_cyclic_variables ( int dim, int **nvr );
/*
 * DESCRIPTION :
 *   Defines the number of variables in each monomial
 *   in the column representation of the cyclic n-roots system.
 *
 * ON ENTRY :
 *   dim     number of equations and total number of variables;
 *   deg     degree of the power series;
 *   nvr     space for dim columns and dim variables in each column.
 *
 * ON RETURN :
 *   nvr     nvr[[i][j] is the number of variables of the j-th monomial
 *           in the i-th column. */

void make_cyclic_columns ( int dim, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Defines the column representation of the cyclic n-roots system.
 *
 * ON ENTRY :
 *   dim     number of equations and total number of variables;
 *   deg     degree of the power series;
 *   nvr     nvr[[i][j] is the number of variables of the j-th monomial
 *           in the i-th column;
 *   idx     space for the variables that appear in every monomial
 *           in every column.
 *
 * ON RETURN :
 *   idx     idx[i][j][k] is the index of the k-th variable which appears
 *           in the j-th monomial of the i-th column. */

void write_cyclic_columns ( int dim, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Writes the exponents of the column representation of
 *   the cyclic n-roots systems. */

int main ( void )
{
   cout << "testing column representation of cyclic n-roots ..." << endl;

   cout << "give the dimension : ";
   int dim; cin >> dim;

   int **nvr = new int*[dim];
   for(int i=0; i<dim; i++) nvr[i] = new int[dim];

   make_cyclic_variables(dim,nvr);

   int ***idx = new int**[dim]; // we have dim columns and
   for(int i=0; i<dim; i++)     // dim monomials in each column
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++) idx[i][j] = new int[nvr[i][j]];
   }
   make_cyclic_columns(dim,nvr,idx);

   write_cyclic_columns(dim,nvr,idx);

   return 0;
}

void make_cyclic_variables ( int dim, int **nvr )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim-1; j++) nvr[i][j] = j+1;
   nvr[0][dim-1] = dim;
   for(int i=1; i<dim; i++) nvr[i][dim-1] = 0;
}

void make_cyclic_columns ( int dim, int **nvr, int ***idx )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=0; k<nvr[i][j]; k++) idx[i][j][k] = (i + k) % dim;
}

void write_cyclic_columns ( int dim, int **nvr, int ***idx )
{
   int *deg = new int[dim];

   for(int i=0; i<dim; i++)
   {
      cout << "column " << i << " :" << endl;
      for(int j=0; j<dim; j++)
      {
         cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";

         for(int k=0; k<dim; k++) deg[k] = 0;
         for(int k=0; k<nvr[i][j]; k++) deg[idx[i][j][k]] = 1;
         for(int k=0; k<dim; k++) cout << " " << deg[k];

         cout << endl;
      }
   }
   free(deg);
}
