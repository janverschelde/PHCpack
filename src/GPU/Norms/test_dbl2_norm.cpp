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
#include "dbl2_norm_kernels.h"
#include "random_vectors.h"
#include "double_double.h"

using namespace std;

void run ( int dim,
           double *vhinorm_device, double *vlonorm_device,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_device, double *wlonorm_device,
           double *whinorm_host, double *wlonorm_host );
/*
  DESCRIPTION :
    Computes norms for random real vectors.

  ON ENTRY :
    dim    dimension of the random vector;

  ON RETURN :
    vhinorm_device  high part of 2-norm compute by the device;
    vlonorm_device  low part of 2-norm computed by the device;
    vhinorm_host    high part of 2-norm compute by the host;
    vlonorm_host    low part of 2-norm computed by the host;
    whinorm_device  high part of 2-norm on device after normalization;
    wlonorm_device  low norm of 2-norm on device after normalization;
    whinorm_host    high part of 2-norm on host after normalization;
    wlonorm_host    low norm of 2-norm on host after normalization.  */

int verify_correctness ( int dim );
/*
 * Computes norms of real vectors of dimension dim
 * and verifies the correctness. */

int main ( void )
{
   int fail = verify_correctness(256);

   return 0;
}

void run ( int dim,
           double *vhinorm_device, double *vlonorm_device,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_device, double *wlonorm_device,
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

   GPU_norm(vhi_device,vlo_device,dim,1,dim,vhinorm_device,vlonorm_device,1);
   GPU_norm(vhi_device,vlo_device,dim,1,dim,whinorm_device,wlonorm_device,1);

   CPU_norm(vhi_host,vlo_host,dim,vhinorm_host,vlonorm_host);
   make_copy(dim,vhi_host,vlo_host,whi_host,wlo_host);
   CPU_normalize(whi_host,wlo_host,dim,*vhinorm_host,*vlonorm_host);
   CPU_norm(whi_host,wlo_host,dim,whinorm_host,vlonorm_host);
}

int verify_correctness ( int dim )
{
   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host
   double vnrm_h[2],vnrm_d[2],wnrm_h[2],wnrm_d[2];

   run(dim,&vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
           &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);

   vnrm_d[0] = vhinorm_device; vnrm_d[1] = vlonorm_device;
   wnrm_d[0] = whinorm_device; wnrm_d[1] = wlonorm_device;
   vnrm_h[0] = vhinorm_host; vnrm_h[1] = vlonorm_host;
   wnrm_h[0] = whinorm_host; wnrm_h[1] = wlonorm_host;

   cout << scientific << setprecision(16);

   cout << "   GPU norm : ";
   dd_write(vnrm_d,32); cout << endl;
   cout << "       after normalization : ";
   dd_write(wnrm_d,32); cout << endl;
   cout << "   CPU norm : ";
   dd_write(vnrm_h,32); cout << endl;
   cout << "       after normalization : ";
   dd_write(wnrm_h,32); cout << endl;

   return 0;
}
