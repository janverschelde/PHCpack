#include <iostream>
#include <iomanip>
#include "mycomplex.h"
#include <qd/qd_real.h>

using namespace std;

int main ( void )
{
   mycomplex<double> a(1.0/3, 1.0/7);
   
   dd_real dd1 = 1.0;
   dd_real dd3 = 3.0;
   dd_real dd7 = 7.0;
   dd_real ddx = dd1/dd3;
   dd_real ddy = dd1/dd7;
   qd_real qd1 = 1.0;
   qd_real qd3 = 3.0;
   qd_real qd7 = 7.0;
   qd_real qdx = qd1/qd3;
   qd_real qdy = qd1/qd7;

   mycomplex<dd_real> b;
   b.real = ddx; b.image = ddy;
   mycomplex<qd_real> c;
   c.real = qdx; c.image = qdy;

   cout << "Testing complex numbers at various precision ...\n" << endl;
   cout << "a : " << a << endl;
   cout << setprecision(32) << "b : " << b << endl;
   cout << setprecision(64) << "c : " << c << endl;

   return 0;
}
