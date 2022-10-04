// The file dbl_bals_flopcounts.cpp defines the functions with prototypes in
// the file dbl_bals_flopcounts.h.

#include <iostream>

void flopcount_dbl_bals_tail
 ( int dim, long long int *add, long long int *mul )
{
   *add += dim*dim; // szt*nbt = dim blocks of threads do dim subtractions
   *mul += dim*dim; // szt*nbt = dim blocks of threads do dim multiplications
}

void flopcount_cmplx_bals_tail
 ( int dim, long long int *add, long long int *mul )
{
   // running szt*nbt dim blocks of threads
   *add += 4*dim*dim; // each block does 4*dim subtractions
   *mul += 4*dim*dim; // each block does 4*dim multiplications
}

void prompt_flopbals_setup
 ( int *seed, int *dim, int *deg, int *szt, int *nbt, int *cdata )
{
   using namespace std;

   cout << "-> give the seed (0 for time) : "; cin >> *seed;
   cout << "-> on complex data ? (0 is no, 1 is yes) : "; cin >> *cdata;
   cout << "-> give the dimension : "; cin >> *dim;
   cout << "-> give the degree of the series : "; cin >> *deg;
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
