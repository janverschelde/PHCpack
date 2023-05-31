// The file cyclic_indices.cpp defines the functions specified in
// the file cyclic_indices.h.

#include <iostream>
#include <cstdlib>
#include "random_polynomials.h"

using namespace std;

void make_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx )
{
   for(int i=0; i<dim-1; i++) nbr[i] = dim;
   nbr[dim-1] = 1;

   for(int i=0; i<dim-1; i++)
   {
      nvr[i] = new int[dim];
      for(int j=0; j<dim; j++) nvr[i][j] = i+1;
   }
   nvr[dim-1] = new int[1]; // 2 monomials in the last polynomial
   nvr[dim-1][0] = dim;     // but constant is stored separately

   for(int i=0; i<dim-1; i++)   // dim monomials in each polynomial
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++)
      {
         idx[i][j] = new int[nvr[i][j]];
         for(int k=0; k<nvr[i][j]; k++) idx[i][j][k] = (j + k) % dim;
         insert_sort(nvr[i][j],idx[i][j]);
      }
   }
   idx[dim-1] = new int*[1];
   idx[dim-1][0] = new int[dim]; // except for the last monomial
   for(int k=0; k<dim; k++) idx[dim-1][0][k] = k;
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
