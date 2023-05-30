/* Tests evaluation and differentiation of the cyclic n-roots system 
 * in double precision with the polynomial system data structure. */

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void write_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Writes the indices of a polynomial system.
 *
 * ON ENTRY :
 *   dim      number of polynomials in the system;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial.    */

int main ( void )
{
   cout << "evaluation and differentiation of cyclic n-roots ..." << endl;

   cout << "-> give the dimension : ";
   int dim; cin >> dim;
   cout << "-> give the degree : ";
   int deg; cin >> deg;

   int *nbr = new int[dim];
   for(int i=0; i<dim-1; i++) nbr[i] = dim;
   nbr[dim-1] = 1;

   int **nvr = new int*[dim]; // number of variables in dim polynomials
   for(int i=0; i<dim-1; i++)
   {
      nvr[i] = new int[dim];
      for(int j=0; j<dim; j++) nvr[i][j] = i+1;
   }
   nvr[dim-1] = new int[1]; // 2 monomials in the last polynomial
   nvr[dim-1][0] = dim;     // but constant is stored separately

   int ***idx = new int**[dim]; // we have dim polynomials and
   for(int i=0; i<dim-1; i++)   // dim monomials in each polynomial
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++)
      {
         idx[i][j] = new int[nvr[i][j]];
         for(int k=0; k<nvr[i][j]; k++) idx[i][j][k] = (j + k) % dim;
      }
   }
   idx[dim-1] = new int*[1];
   idx[dim-1][0] = new int[dim]; // except for the last monomial
   for(int k=0; k<dim; k++) idx[dim-1][0][k] = 1;

   write_polynomial_indices(dim,nbr,nvr,idx);

   return 0;
}

void write_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx )
{
   for(int i=0; i<dim; i++)
   {
      cout << "polynomial " << i
           << " has " << nbr[i] << " monomials :" << endl;
      for(int j=0; j<nbr[i]; j++)
      {
         for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
         cout << endl;
      }
   }
}
