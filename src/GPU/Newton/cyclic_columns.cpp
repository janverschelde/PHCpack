// The file cyclic_columns.cpp defines the functions specified in
// the file cyclic_columns.h

#include <iostream>
#include <cstdlib>
#include "cyclic_columns.h"

using namespace std;

void make_cyclic_variables ( int nbrcol, int dim, int **nvr )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim-1; j++) nvr[i][j] = j+1;
   nvr[0][dim-1] = dim; // only first column has nonempty monomial
   for(int i=1; i<nbrcol; i++) nvr[i][dim-1] = 0;
}

void make_cyclic_columns ( int nbrcol, int dim, int **nvr, int ***idx )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         for(int k=0; k<nvr[i][j]; k++) idx[i][j][k] = (i + k) % dim;
}

void write_cyclic_columns ( int nbrcol, int dim, int **nvr, int ***idx )
{
   int *deg = new int[dim];

   for(int i=0; i<nbrcol; i++)
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
