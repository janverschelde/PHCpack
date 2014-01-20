#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <assert.h>
#include <cmath>
#include "vector_functions.h"
#include "gqd_type.h"
#include "gqd_qd_utilT.h"
#include "norm_kernels.h"
#include "norm_host.h"
#include "DefineType.h"

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

   // vector on the host is stored in v
   complexH<T1>*v = new complexH<T1>[dim];
   // v_h contains the GPU complex vector, as stored on the host
   complex<T>* v_h = new complex<T>[dim];
   random_point(dim,v_h,v);  // random numbers on unit circle

   // for normalization, we make copies of v and v_h in w and w_h :
   complexH<T1>*w = new complexH<T1>[dim];

   T vnorm_h,wnorm_h;
   T1 vnorm,wnorm;

   // GPU computation of the norm
   if(mode==0 || mode==2)
   {
      GPU_norm(v_h,dim,1,BS,&vnorm_h);
      GPU_norm(v_h,dim,freq,BS,&wnorm_h);
   }

   // CPU computation of the norm
   if(mode==1 || mode==2)
      for(int i=0; i<freq; i++)
      {
         CPU_norm(v,dim,&vnorm);
         make_copyH(dim,v,w);
         CPU_normalize(w,dim,vnorm);
         CPU_norm(w,dim,&wnorm);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      T1 temp;
      gqd2qd(&vnorm_h,&temp);
      cout << "sqrt(" << dim << ") = " << sqrt(dim) << endl;
      cout << "GPU norm : " << temp << endl;
      gqd2qd(&wnorm_h,&temp);
      cout << "GPU norm after normalization : " << temp << endl;
      cout << "CPU norm : " << vnorm << endl;
      cout << "CPU norm after normalization : " << wnorm << endl;
   }

   return 0;
}
