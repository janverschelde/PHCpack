/* The file prompt_for_setup.cpp contains the definition of a function. */

#include <iostream>

using namespace std;

void prompt_baqr_setup
 ( int *seed, int *szt, int *nbt, int *nrows, int *vrb, int *mode )
{
   cout << "Give the seed (0 for time) : "; cin >> *seed;

   cout << "Give the size of each tile : "; cin >> *szt;

   cout << "Give the number of tiles : "; cin >> *nbt;

   const int ncols = (*szt)*(*nbt);
   bool done = false;

   while(!done)
   {
      cout << "Give the number of rows (" << " >= " << ncols << " ) : ";
      cin >> *nrows;
      if(*nrows >= ncols)
         done = true;
      else
      {
         cout << "Number of rows = " << *nrows
              << " < " << ncols
              << " = number of columns." << endl;
         cout << "Please try again." << endl;
      }
   }
   cout << "Give the verbose level : "; cin >> *vrb;

   cout << "Enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> *mode;
}
