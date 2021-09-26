/* The file dbl_data_files.cpp defines the functions specified in
 * the file dbl_data_files.h. */

#include <iostream>
#include <iomanip>
#include <fstream>
#include "dbl_data_files.h"

void dbl_write_matrix ( string name, int dim, double **A )
{
   ofstream outs(name.c_str());

   outs << scientific << setprecision(16);

   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++) outs << " " << A[i][j];
      outs << endl;
   }
   outs.close();
}

void cmplx_write_matrix ( string name, int dim, double **Are, double **Aim )
{
   ofstream outs(name.c_str());

   outs << scientific << setprecision(16);

   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++)
         outs << " " << Are[i][j] << " " << Aim[i][j];

      outs << endl;
   }
   outs.close();
}

void dbl_read_matrix ( string name, int dim, double **A )
{
   ifstream infs(name.c_str());

   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++) infs >> A[i][j];
   }
   infs.close();
}

void cmplx_read_matrix ( string name, int dim, double **Are, double **Aim )
{
   ifstream infs(name.c_str());

   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++)
      {
         infs >> Are[i][j];
         infs >> Aim[i][j];
      }
   }
   infs.close();
}
