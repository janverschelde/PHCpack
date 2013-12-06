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

int parse_arguments
 ( int argc, char *argv[],
   int *BS, int *dim, int *r, int *mode )
{
   if(argc < 5)
   {
      cout << argv[0] << "needs four parameters" << endl;
      cout << "please try again..." << endl; return 1;
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

