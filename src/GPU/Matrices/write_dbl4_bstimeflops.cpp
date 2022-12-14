/* Defines the function with prototype in write_dbl4_bstimeflops. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "write_dbl4_bstimeflops.h"

void write_dbl4_bstimeflops
 ( int sizetile, int numtiles, int ctype,
   double invlapsed, double mullapsed, double sublapsed, double elapsedms,
   double timelapsed, long int addcnt, double addover,
   long int mulcnt, double mulover, long int divcnt, double divover )
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
   cout << "             Number of additions/subtractions : ";
   if(addover == 0.0)
      cout << addcnt << " x 89 " << endl;
   else
   {
      cout << scientific << setprecision(16);
      addover = addover + (double) addcnt;
      cout << addover << " x 89 " << endl;
   }
   cout << "                    Number of multiplications : ";
   if(mulover == 0.0)
      cout << mulcnt << " x 336 " << endl;
   else
   {
      cout << scientific << setprecision(16);
      mulover = mulover + (double) mulcnt;
      cout << mulover << " x 336 " << endl;
   }
   cout << "                          Number of divisions : ";
   if(divover == 0.0)
      cout << divcnt << " x 893 " << endl;
   else
   {
      cout << scientific << setprecision(16);
      divover = divover + (double) divcnt;
      cout << divover << " x 893 " << endl;
   }
   double kernflops,wallflops,theflopcnt;

   cout << "    Total number of floating-point operations : ";
   if((addover == 0.0) && (mulover == 0.0) && (divover == 0.0))
   {
      long long int flopcnt = 89*addcnt + 336*mulcnt + 893*divcnt;
      cout << flopcnt << endl;
      cout << endl;
      kernflops = 1000.0*((double) flopcnt)/elapsedms;
      wallflops = ((double) flopcnt)/timelapsed;
      theflopcnt = (double) flopcnt;
   }
   else
   {
      double flopcnt = 89*addover + 336*mulover + 893*divover;    
      if(addover == 0.0) flopcnt += 89*addcnt;
      if(mulover == 0.0) flopcnt += 336*mulcnt;
      if(divover == 0.0) flopcnt += 893*divcnt;
      cout << flopcnt << endl;
      cout << endl;
      kernflops = 1000.0*flopcnt/elapsedms;
      wallflops = flopcnt/timelapsed;
      theflopcnt = flopcnt;
   }
   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 4*sizetile*numtiles*(numtiles+1)*4 + 8*sizetile*numtiles*4;
   else
      bytecnt = 4*sizetile*numtiles*(numtiles+1)*8 + 8*sizetile*numtiles*8;

   cout << "    Total number of bytes : " << bytecnt << endl << endl;

   const int gigacnt = pow(2.0,30);
   double intensity = theflopcnt/bytecnt;

   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl << endl;

   cout << "  Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << "   Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
}
