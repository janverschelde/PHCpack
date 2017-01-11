#include <iostream>
#include <iomanip>
#include "mycomplex.h"
#include <qd/qd_real.h>

using namespace std;

template <class real>
mycomplex<real> square_root ( mycomplex<real> x, int prc )
{
   mycomplex<real> wrk,val,dx;
   cout << setprecision(prc) << "\nx : " << x;
   wrk = x;
   for(int i=0; i<8; i++)
   {
      val = wrk*wrk - x;
      dx = val/(wrk+wrk);
      cout << setprecision(4) << "dx : " << dx;
      wrk = wrk - dx;
   }
   return wrk;
}

int main ( void )
{
   char choice;

   cout << "Testing complex numbers at various precision ...\n" << endl;
   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   if(choice == '0')
   {
      mycomplex<double> a(1.0/3, 1.0/7);
      mycomplex<double> r;
      r = square_root<double>(a, 16);
      cout << setprecision(16) << "\nr : " << r;
      cout << setprecision(16) << "r*r : " << r*r;
   }
   else if(choice == '1')
   {
      dd_real dd1 = 1.0;
      dd_real dd3 = 3.0;
      dd_real dd7 = 7.0;
      dd_real ddx = dd1/dd3;
      dd_real ddy = dd1/dd7;
      mycomplex<dd_real> b,r;
      b.re = ddx; b.im = ddy;
      r = square_root<dd_real>(b, 32);
      cout << setprecision(32) << "\nr : " << r;
      cout << setprecision(32) << "r*r : " << r*r;
   }
   else if(choice == '2')
   {
      qd_real qd1 = 1.0;
      qd_real qd3 = 3.0;
      qd_real qd7 = 7.0;
      qd_real qdx = qd1/qd3;
      qd_real qdy = qd1/qd7;
      mycomplex<qd_real> c,r;
      c.re = qdx; c.im = qdy;
      r = square_root<qd_real>(c, 64);
      cout << setprecision(64) << "\nr : " << r;
      cout << setprecision(64) << "r*r : " << r*r;
   }
   else
      cout << "Invalid choice for the precision." << endl;

   return 0;
}
