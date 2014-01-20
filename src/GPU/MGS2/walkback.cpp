#include <iostream>

using namespace std;

int main ( void )
{
   cout << "give the dimension : ";
   int dim; cin >> dim;
   cout << "give the block size : ";
   int BS; cin >> BS;

   if(BS == dim)
   {
      for(int k=dim-1; k>=0; k--) 
      {
         int j = 0;
         int ind = (dim - k)*(dim + 3 + k)/2 + j; 
         cout << "k = " << k << "  ind = " << ind << endl;
         for(j=0; j<BS; j++)
            if(j <= k)
            {
               cout << "thread " << j << " :";
               cout << "  ind + j = " << ind + j << endl;
            }
      }
   }
   else
   {
      int rnd = dim/BS;
      cout << "number of rounds : " << rnd << endl;
      for(int p = rnd-1; p>=0; p--)
      {
         int offset = p*BS;
         cout << "round " << p << " offset = " << offset << endl;
         cout << "*** indexing for block 0 ***" << endl;
         for(int k=BS-1; k>=0; k--) 
         {
            int j = 0;
            int ind = (dim - k-offset)*(dim + 3 + k+offset)/2 + j; 
            cout << "k = " << k << "  ind + offset = " << ind+offset << endl;
            for(j=0; j<BS; j++)
               if(j <= k) 
               {
                  cout << "thread " << j << " :";
                  cout << "  ind + offset + j = " << ind + offset + j << endl;
               }
         }
         for(int b=1; b<=p; b++) 
         {
            cout << "*** indexing for block " << b << " ***" << endl;
            int block_offset = b*BS;
            cout << "block offset = " << block_offset << endl;
            for(int k=BS-1; k>=0; k--) 
            {
               int j = 0;
               int ind = (dim-k-offset)*(dim+3+k+offset)/2 + j; 
               cout << "k = " << k << "  ind + offset - block_offset = ";
               cout << ind + offset - block_offset << endl;
               for(j=0; j<BS; j++)
               {
                  cout << "thread " << j << " :";
                  cout << "  ind + offset - block_offset + j = ";
                  cout << ind + offset - block_offset + j << endl;
               }
            }
         }
      }
   }
}
