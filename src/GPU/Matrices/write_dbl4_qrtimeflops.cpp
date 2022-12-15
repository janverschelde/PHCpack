/* Defines the function with prototype in write_dbl4_qrtimeflops. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "write_dbl4_qrtimeflops.h"

void write_dbl4_qrtimeflops
 ( int ctype, int nrows, int ncols,
   double houselapsedms, double RTvlapsedms, double tileRlapsedms,
   double vb2Wlapsedms, double WYTlapsedms, double QWYTlapsedms,
   double Qaddlapsedms, double YWTlapsedms, double YWTClapsedms,
   double Raddlapsedms, double timelapsed,
   long long int addcnt, long long int mulcnt,
   long long int divcnt, long long int sqrtcnt )
{
   using namespace std;

   cout << fixed << setprecision(3);

   cout << "         Time spent by the Householder kernel : "
        << houselapsedms << " milliseconds." << endl;
   cout << "      Time spent by the kernel for beta*R^T*v : "
        << RTvlapsedms << " milliseconds." << endl;
   cout << "  Time spent by the kernel to reduce one tile : "
        << tileRlapsedms << " milliseconds." << endl;
   cout << "    Time spent by the kernel for the W matrix : "
        << vb2Wlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing Y*W^T : "
        << YWTlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing Q*WYT : "
        << QWYTlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing YWT*C : "
        << YWTClapsedms << " milliseconds." << endl;
   cout << "Time spent by the kernel for adding QWYT to Q : "
        << Qaddlapsedms << " milliseconds." << endl;
   cout << "Time spent by the kernel for adding R to YWTC : "
        << Raddlapsedms << " milliseconds." << endl;

   const double totlapsedms = houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + YWTlapsedms + QWYTlapsedms
      + YWTClapsedms + Qaddlapsedms + Raddlapsedms;

   cout << "                    Time spent by all kernels : "
        << totlapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed << " seconds." << endl;
   cout << endl;
   cout << "             Number of additions/subtractions : "
        << addcnt << " x 89 " << endl;
   cout << "                    Number of multiplications : "
        << mulcnt << " x 336 " << endl;
   cout << "                          Number of divisions : "
        << divcnt << " x 893 " << endl;
   cout << "                    Number of calls to sqrt() : "
        << sqrtcnt << " x 1345 " << endl;

   long long int flopcnt = 89*addcnt + 336*mulcnt
                         + 893*divcnt + 1345*sqrtcnt;

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;
   cout << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 4*nrows*ncols + 4*nrows*nrows;
   else
      bytecnt = 8*nrows*ncols + 8*nrows*nrows;

   cout << "    Total number of bytes : " << bytecnt << endl << endl;

   double intensity = ((double) flopcnt)/bytecnt;

   cout << "     Arithmetic intensity : "
        << scientific << setprecision(5) << intensity
        << " #flops/#bytes" << endl << endl;

   double kernflops = 1000.0*((double) flopcnt)/totlapsedms;
   double wallflops = ((double) flopcnt)/timelapsed;
   const int gigacnt = pow(2.0,30);

   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
}
