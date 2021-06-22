/* Tests the jobs to solve an upper triangular system. */

#include <iostream>

using namespace std;

void list_upper_jobs ( int dim );
/*
 * Writes all jobs to solve an upper triangular system
 * of dimension dim. */

int main ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Jobs to solve a system of dimension "
        << dim << " : " << endl;

   list_upper_jobs(dim);

   return 0;
}

void list_upper_jobs ( int dim )
{
   int lvl;
   const int maxlvl = 3*dim-2;
   int* freqops = new int[maxlvl];

   for(int i=0; i<maxlvl; i++) freqops[i] = 0;

   cout << "Level of operations in the backward substitution :" << endl;

   for(int i=dim-1; i>=0; i--)
   {
      cout << "prod = b[" << i << "], level 0" << endl;

      lvl = 1;

      for(int j=dim-1; j>i; j--)
      {
         cout << "work = U[" << i << "][" << j << "]*x[" << j
              << "], level " << ++lvl << endl;
         freqops[lvl] = freqops[lvl] + 1;

         cout << "prod = prod - work, level = " << ++lvl << endl;
         freqops[lvl] = freqops[lvl] + 1;

         lvl = lvl + 1;
      }
      cout << "work = 1/U[" << i << "][" << i << "], level 0" << endl;
      freqops[0] = freqops[0] + 1;

      cout << "x[" << i << "] = prod/U["
           << i << "][" << i << "], level " << lvl << endl;
      freqops[lvl] = freqops[lvl] + 1;
   }
   cout << "The frequency table : " << endl;
   for(int i=0; i<maxlvl; i++)
      cout << i << " : " << freqops[i] << endl;
}
