// The file prompt_newton_setup.cpp defines the function with prototype in
// the file prompt_newton_setup.h.

#include <iostream>
#include <time.h>
#include "unimodular_matrices.h"

using namespace std;

void prompt_newton_setup
 ( int *seed, int *szt, int*nbt, int *dim, int *deg, int *size, int *vrblvl,
   int *mode, int *nbritr, int *nbrcol, int *nbsteps, int *cdata )
{
   cout << "-> give the seed (0 for time) : "; cin >> *seed;
   cout << "-> on complex data ? (0 is no, 1 is yes) : "; cin >> *cdata;

   int posvals = 1;

   prompt_dimensions(dim,deg,size,&posvals,vrblvl,nbritr,nbsteps);

   *nbrcol = 1;

   if(*nbritr == -4) *nbrcol = 2; // 2-column lower/upper triangle

   if(*nbritr == -3)
   {
      *nbrcol = 0;
      while(*nbrcol < 1)
      {
         cout << "-> give the number of columns : "; cin >> *nbrcol;

         if(*nbrcol < 1)
            cout << "-> number of columns must be at least 1, retry."
                 << endl;

         if(*nbrcol > *dim)
         {
            cout << "-> number of columns must not be more than "
                 << *dim << ", retry." << endl;
            *nbrcol = 0;
         }
      }
   }
   cout << "-> enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> *mode;

   if(*mode != 1)
   {
      cout << "-> give the number of tiles : "; cin >> *nbt;
      cout << "-> give the size of each tile : "; cin >> *szt;
      int p = (*szt)*(*nbt);

      while(p != *dim)
      {
          cout << "Dimension = " << *dim << " != " << *szt << " * " << *nbt
               << ", retry." << endl;
          cout << "-> give the size of each tile : "; cin >> *szt;
          cout << "-> give the number of tiles : "; cin >> *nbt;
          p = (*szt)*(*nbt);
      }
   }
}
