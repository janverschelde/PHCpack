/* Defines the function with prototype in write_dbl8_bstimeflops. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "write_dbl8_bstimeflops.h"

void write_dbl8_bstimeflops
 ( int sizetile, int numtiles, int ctype,
   double invlapsed, double mullapsed, double sublapsed, double elapsedms,
   double timelapsed, long int addcnt, long int mulcnt, long int divcnt )
{
   using namespace std;

   cout << fixed << setprecision(3);

   cout << "          Time spent to invert diagonal tiles : ";
   cout << invlapsed << " milliseconds." << endl;
   cout << "   Time spent to multiply with inverted tiles : ";
   cout << mullapsed << " milliseconds." << endl;
   cout << "             Time spent for back substitution : ";
   cout << sublapsed << " milliseconds." << endl;
   cout << "                    Time spent by all kernels : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed << " seconds." << endl;
   cout << endl;
   cout << "             Number of additions/subtractions : "
        << addcnt << " x 270 " << endl;
   cout << "                    Number of multiplications : "
        << mulcnt << " x 1742 " << endl;
   cout << "                          Number of divisions : "
        << divcnt << " x 5126 " << endl;

   long long int flopcnt = 270*addcnt + 1742*mulcnt + 5126*divcnt;
   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;
   cout << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 4*sizetile*numtiles*(numtiles+1)*8 + 8*sizetile*numtiles*8;
   else
      bytecnt = 4*sizetile*numtiles*(numtiles+1)*16 + 8*sizetile*numtiles*16;

   cout << "    Total number of bytes : " << bytecnt << endl << endl;

   const int gigacnt = pow(2.0,30);
   double intensity = ((double) flopcnt)/bytecnt;
   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl << endl;

   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
   double wallflops = ((double) flopcnt)/timelapsed;

   cout << "  Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << "   Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
}
