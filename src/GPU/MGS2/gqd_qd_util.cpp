// Conversions between real multiprecision CPU and GPU types
// and parsing command line arguments as defined in gqd_qd_util.h.

#include "gqd_type.h"
#include "vector_types.h"
#include <qd/qd_real.h>
#include <cstdlib>
#include <iostream>

using namespace std;

void qd2gqd ( qd_real *a, gqd_real *b )
{
   b->x = a->x[0];
   b->y = a->x[1];
   b->z = a->x[2];
   b->w = a->x[3];
}

void gqd2qd ( gqd_real *a, qd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
   b->x[2] = a->z;
   b->x[3] = a->w;
}

void qd2gqd ( dd_real *a, gdd_real *b )
{
   b->x = a->x[0];
   b->y = a->x[1];
}

void gqd2qd ( gdd_real *a, dd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
}

void qd2gqd ( double *a, double *b )
{
   *b = *a;
}

void gqd2qd ( double *a, double *b )
{
   *b = *a;
}

int print_modes ( void )
{
   cout << "The modes to the modified Gram-Schmidt method are\n";
   cout << " 0 : GPU accelerated orthonormalization of complex vectors\n";
   cout << " 1 : CPU computed orthonormalization of complex vectors\n";
   cout << " 2 : GPU+CPU complex orthonormalization with checks\n";
   cout << " 3 : GPU accelerated QR decomposition of complex vectors\n";
   cout << " 4 : CPU computed QR decomposition of complex vectors\n";
   cout << " 5 : GPU+CPU complex QR decomposition with checks\n";
   cout << " 6 : GPU accelerated complex least squares solving\n";
   cout << " 7 : CPU computed complex least squares solving\n";
   cout << " 8 : GPU+CPU complex least squares solving with checks\n";
   cout << " 9 : GPU accelerated complex Newton on H-Chandrasekhar\n";
   cout << "10 : CPU computed complex Newton on H-Chandrasekhar\n";
   cout << "11 : GPU+CPU complex Newton on H-Chandrasekhar with checks\n";
   cout << "12 : GPU accelerated orthonormalization of real vectors\n";
   cout << "13 : CPU computed orthonormalization of real vectors\n";
   cout << "14 : GPU+CPU orthonormalization of real vectors with checks\n";
   cout << "15 : GPU accelerated QR decomposition of real vectors\n";
   cout << "16 : CPU computed QR decomposition of real vectors\n";
   cout << "17 : GPU+CPU QR decomposition of real vectors with checks\n";
   cout << "18 : GPU accelerated real least squares solving\n";
   cout << "19 : CPU computed real least squares solving\n";
   cout << "20 : GPU+CPU real least squares solving with checks\n";
   cout << "21 : GPU accelerated real Newton on H-Chandrasekhar\n";
   cout << "22 : CPU computed real Newton on H-Chandrasekhar\n";
   cout << "23 : GPU+CPU real Newton on H-Chandrasekhar with checks\n";
   cout << "24 : GPU full complex Newton on H-Chandrasekhar\n";
   cout << "25 : GPU full complex Newton on H-Chandrasekhar with output\n";
}

int parse_arguments
 ( int argc, char *argv[],
   int *BS, int *dim, int *r, int *mode )
{
   if(argc < 5)
   {
      cout << argv[0] << " needs four parameters, for example:" << endl;
      cout << argv[0] << " BS dim r mode" << endl;
      cout << " where BS is the number of threads in a block," << endl;
      cout << "       dim is the dimension of the problem," << endl;
      cout << "       r is the number of repeated runs, and" << endl;
      cout << "       mode is the execution mode." << endl;
      print_modes();
      cout << "Please try again ..." << endl; return 1;
   }
   *BS = atoi(argv[1]);     // block size
   *dim = atoi(argv[2]);    // dimension
   *r = atoi(argv[3]);      // number of repeated runs
   *mode = atoi(argv[4]);   // execution mode
   return 0;
}

double random_double ( void )
{
   return rand()/((double) RAND_MAX);
}

