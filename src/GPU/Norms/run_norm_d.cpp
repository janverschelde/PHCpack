// This is the main program.

// Including the vector_types.h is needed for the include of
// the headers for the dbl_norm_kernels.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_norm_kernels.h"
#include "dbl_norm_host.h"

using namespace std;

int parse_arguments
 ( int argc, char *argv[], int *blocksize, int *dim, int *freq, int *mode );
/*
   Parses the argc arguments on the command line in argv[].
   Returns 0 if four numbers are given, 1 otherwise for failure.

   ON RETURN :
     blocksize  the block size, number of threads in a block;
     dim        dimension of the vectors;
     freq       frequency of the runs;
     mode       execution mode is 0, 1, or 2
                0 : GPU run, no output,
                1 : CPU run, no output,
                2 : both GPU and CPU run with output. */

double random_double ( void );
// Returns a random double in [-1,+1].

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

   double *v_host = new double[dim];   // vector on the host
   double* v_device = new double[dim]; // vector on the device
   for(int k=0; k<dim; k++)
   {
      double r = random_double();
      v_host[k] = r;
      v_device[k] = r;
   }

   double *w_host = new double[dim];   // a copy for normalization

   double vnorm_device,wnorm_device;
   double vnorm_host,wnorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(v_device,dim,1,BS,&vnorm_device);
      GPU_norm(v_device,dim,freq,BS,&wnorm_device);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(v_host,dim,&vnorm_host);
         make_copy(dim,v_host,w_host);
         CPU_normalize(w_host,dim,vnorm_host);
         CPU_norm(w_host,dim,&wnorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << "GPU norm : " << vnorm_device << endl;
      cout << "GPU norm after normalization : " << wnorm_device << endl;
      cout << "CPU norm : " << vnorm_host << endl;
      cout << "CPU norm after normalization : " << wnorm_host << endl;
   }

   return 0;
}

int parse_arguments
 ( int argc, char *argv[], int *blocksize, int *dim, int *freq, int *mode )
{
   if(argc < 5)
   {
      cout << argv[0] << " needs four parameters, for example:" << endl;
      cout << argv[0] << " blocksize dim freq mode" << endl;
      cout << " where blocksize is the number of threads in a block," << endl;
      cout << "       dim is the dimension of the problem," << endl;
      cout << "       freq is the number of repeated runs, and" << endl;
      cout << "       mode is 0, 1, or 2 for execution mode." << endl;
      cout << "Please try again ..." << endl;

      return 1;
   }
   else
   {
      *blocksize = atoi(argv[1]);  // number of threads in a block
      *dim = atoi(argv[2]);        // dimension
      *freq = atoi(argv[3]);       // number of repeated runs
      *mode = atoi(argv[4]);       // execution mode

      return 0;
   }
}

double random_double ( void )
{
   double r = (double) rand();

   r = r/RAND_MAX; // r is in [0,1]

   return (2.0*r - 1.0); // result is in [-1, 1]
}
