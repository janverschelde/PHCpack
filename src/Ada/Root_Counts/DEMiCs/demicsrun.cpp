#include <iostream>
#include "inputData.h"
#include "demicsrun.h"

using namespace std;

extern "C" int demicsrun
 ( int verbose, int dimension, int nsupports,
   int* mixtype, int* cardsup, int* coordinates )
{
   if(verbose == 1)
      write_data(dimension,nsupports,mixtype,cardsup,coordinates);

   dataSet Data;
   Data.Dim = dimension;
   Data.supN = nsupports;
   Data.type = new int[nsupports];
   for(int k=0; k<Data.supN; k++) Data.type[k] = mixtype[k];
   Data.termSet = new int[nsupports];
   for(int k=0; k<Data.supN; k++) Data.termSet[k] = cardsup[k];

   if(verbose == 1)
   {
      cout << endl;
      cout << "The dimension, the number of distinct support sets," << endl;
      cout << "the number of points in each support set, and" << endl;
      cout << "the number of occurrences of each support set :" << endl;
      cout << endl;

      Data.info_preamble();
   }

   Data.termSumNum = 0;
   for(int k=0; k<Data.supN; k++)
      Data.termSumNum = Data.termSumNum + Data.termSet[k];

   Data.support = new double[Data.termSumNum*Data.Dim];

   int offset = 0;
   int idxcoordinates = 0;
   for(int i=0; i<Data.supN; i++)
   {
      for(int j=0; j<Data.termSet[i]; j++)
      { 
         for(int k=0; k<Data.Dim; k++)
         {
            double x = coordinates[idxcoordinates++];
            Data.support_in(offset+j,k,x);
         }
      }
      offset = offset + Data.termSet[i];
   }

   if(verbose == 1)
   {
      cout << endl;
      cout << "The points in the support sets : " << endl << endl;

      Data.info_supports();
   }
   return 0;
}

void write_data
 ( int dimension, int nsupports,
   int* mixtype, int* cardsup, int* coordinates )
{
   cout << "The dimension : " << dimension << endl;
   cout << "The number of supports : " << nsupports << endl;

   cout << "Mixture type :";
   for(int k=0; k<nsupports; k++) cout << " " << mixtype[k];
   cout << endl;

   cout << "Cardinalities :";
   for(int k=0; k<nsupports; k++) cout << " " << cardsup[k];
   cout << endl;

   cout << "The points in the support sets ";
   int idx = 0;
   for(int i=0; i<nsupports; i++)
   {
      cout << endl;
      for(int j=0; j<cardsup[i]; j++) 
      {
         for(int k=0; k<dimension; k++)
         {
            cout << " " << coordinates[idx];
            idx = idx + 1;
         }
         cout << endl;
      }
   }
}
