/* Defines the function with prototype in write_dbl_bstimeflops. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "write_dbl_bstimeflops.h"

void write_dbl_bstimeflops
 ( double invlapsed, double mullapsed, double sublapsed, double elapsedms,
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
        << addcnt << endl;
   cout << "                    Number of multiplications : "
        << mulcnt << endl;
   cout << "                          Number of divisions : "
        << divcnt << endl;

   long long int flopcnt = addcnt + mulcnt + divcnt;

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;
   cout << endl;
   cout << scientific << setprecision(3);

   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
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
