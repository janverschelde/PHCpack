#include <iostream>
#include "global.h"
#include "inputData.h"
#include "mvc.h"
#include "demicsrun.h"

using namespace std;

extern "C" int demicsrun
 ( int verbose, int dimension, int nsupports,
   int* mixtype, int* cardsup, int* coordinates )
{
   if(verbose == 1)
      write_data(dimension,nsupports,mixtype,cardsup,coordinates);

   dataSet Data;

   fill_preamble(Data,verbose,dimension,nsupports,mixtype,cardsup);
   fill_supports(Data,verbose,coordinates);
   fill_complete(Data,verbose);

   mvc MV_Comp;

   MV_Comp.allocateAndIni(Data,1,verbose);
   MV_Comp.Enum();

   return 0;
}

extern "C" int demicsfly
 ( int verbose, int dimension, int nsupports,
   int* mixtype, int* cardsup, int* coordinates, double* lifvals )
{
   if(verbose == 1)
      write_fly_data(dimension,nsupports,mixtype,cardsup,coordinates,lifvals);

   dataSet Data;

   fill_preamble(Data,verbose,dimension,nsupports,mixtype,cardsup);
   fill_supports(Data,verbose,coordinates);
   fill_complete(Data,verbose);

   mvc MV_Comp;

   MV_Comp.initialize_with_lifting(Data,lifvals,1,verbose);
   MV_Comp.Enum();

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

void write_fly_data
 ( int dimension, int nsupports,
   int* mixtype, int* cardsup, int* coordinates, double* lifvals )
{
   write_data(dimension,nsupports,mixtype,cardsup,coordinates);

   cout << "The lifting values for the points in the support sets ";
   int idx = 0;
   for(int i=0; i<nsupports; i++)
   {
      cout << endl;
      for(int j=0; j<cardsup[i]; j++) 
      {
         cout << " " << lifvals[idx];
         idx = idx + 1;
      }
   }
}

void fill_preamble
 ( dataSet& Data,
   int verbose, int dimension, int nsupports, int* mixtype, int* cardsup )
{
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
}

void fill_supports ( dataSet& Data, int verbose, int* coordinates )
{
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
}

void fill_complete ( dataSet& Data, int verbose )
{
   Data.typeMax = Data.type[0];
   for(int k=1; k<Data.supN; k++)
      if(Data.type[k] > Data.typeMax) Data.typeMax = Data.type[k];

   if(verbose == 1) cout << "typeMax : " << Data.typeMax << endl;

   Data.termMax = Data.termSet[0];
   for(int k=1; k<Data.supN; k++)
      if(Data.termSet[k] > Data.termMax) Data.termMax = Data.termSet[k];

   if(verbose == 1) cout << "termMax : " << Data.termMax << endl;

   Data.termStart = new int[Data.supN+1];
   Data.termStart[0] = 0;
   int totalsum = 0;
   for(int k=1; k<Data.supN+1; k++)
   {
      totalsum += Data.termSet[k-1];
      Data.termStart[k] = totalsum;
   }

   if(verbose == 1)
   {
      cout << "termStart =";
      for(int k=0; k<Data.supN+1; k++) cout << " " << Data.termStart[k];
      cout << endl;
      cout << "termSumNum : " << Data.termSumNum << endl;
   }
   Data.coef = new double[2*Data.termSumNum];
}
