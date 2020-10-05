/* This program runs automatic tests to compute norms of real vectors
   in double double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl2_norm_host.h"
#include "random_vectors.h"

using namespace std;

void run ( int dim,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_host, double *wlonorm_host );
/*
  DESCRIPTION :
    Computes norms for random real vectors.

  ON ENTRY :
    dim    dimension of the random vector;

  ON RETURN :
    vhinorm_host    high part of 2-norm compute by the host;
    vlonorm_host    low part of 2-norm computed by the host;
    whinorm_host    high part of 2-norm on host after normalization;
    wlonorm_host    low norm of 2-norm on host after normalization.  */

int verify_correctness ( int dim );
/*
 * Computes norms of real vectors of dimension dim
 * and verifies the correctness. */

int main ( void )
{
   int fail = verify_correctness(64);

   return 0;
}

void run ( int dim,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_host, double *wlonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vhi_host = new double[dim];   // high parts on the host
   double* vlo_host = new double[dim];   // low parts on the host
   double* vhi_device = new double[dim]; // high parts on the device
   double* vlo_device = new double[dim]; // low parts on the device
   double* whi_host = new double[dim];   // high parts copy for normalization
   double* wlo_host = new double[dim];   // low parts copy for normalization

   random_double2_vectors(dim,vhi_host,vlo_host,vhi_device,vlo_device);

   CPU_norm(vhi_host,vlo_host,dim,vhinorm_host,vlonorm_host);
   make_copy(dim,vhi_host,vlo_host,whi_host,wlo_host);
   CPU_normalize(whi_host,wlo_host,dim,*vhinorm_host,*vlonorm_host);
   CPU_norm(whi_host,wlo_host,dim,whinorm_host,vlonorm_host);
}

int verify_correctness ( int dim )
{
   double vhinorm_host,vlonorm_host; // norm before normalization
   double whinorm_host,wlonorm_host; // norm after normalization

   run(dim,&vhinorm_host,&vlonorm_host,&whinorm_host,&wlonorm_host);

   cout << scientific << setprecision(16);
   cout << "   CPU norm : " << vhinorm_host << endl;
   cout << "       after normalization : " << whinorm_host << endl;
}
