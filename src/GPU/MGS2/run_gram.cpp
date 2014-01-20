// This is the main program to compute the Gram matrix 
// of a sequence of random vectors of complex numbers on the unit circle.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "vector_functions.h"
#include "gqd_type.h"
#include "gqd_qd_util.h"
#include "DefineType.h"
#include "gram_host.h"
#include "gram_kernels.h"

using namespace std;

int main ( int argc, char *argv[] )
{
   // initialization of the execution parameters

   int BS,dim,freq,mode;
   if(parse_arguments(argc,argv,&BS,&dim,&freq,&mode) == 1) return 1;

   int timevalue;
   if(mode == 2)
      timevalue = time(NULL); // no fixed seed to verify correctness
   else
      timevalue = 1287178355; // fixed seed for timings
   srand(timevalue);

   // the sequence of vectors on the host is stored in v
   complexH<T1>** v = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++) v[i] = new complexH<T1>[dim];
   // v_h is for GPU processing, stored on the host, as one long array
   complex<T>* v_h = new complex<T>[dim*dim];
   
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         double temp = random_double()*2*M_PI;
         v[i][j].init(cos(temp),sin(temp));
         v_h[i*dim+j].initH(cos(temp),sin(temp));
      }
   // allocating memory for the Gram matrix g and g_h
   complexH<T1>** g = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++) g[i] = new complexH<T1>[dim];
   complex<T>* g_h = new complex<T>[dim*dim];

   if(mode == 0 || mode == 2) GPU_gram(v_h,g_h,dim,freq,BS);

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<freq; i++) CPU_gram(v,g,dim);
      // if(mode == 2) print_gram_matrices(g,g_h,dim); // only for small dim
      if(mode == 2) print_difference(g,g_h,dim);
   }

   return 0;
}
