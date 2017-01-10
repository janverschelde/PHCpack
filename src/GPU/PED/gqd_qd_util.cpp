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
   int *BS, int *dim, int *NM, int *NV, int *deg, int *r, int *mode )
/* Parses the arguments on the command line.
   Returns 0 if okay, otherwise the function returns 1. */
{
   if(argc < 8)
   {
      cout << argv[0] << " needs seven parameters :" << endl;
      cout << "(1) block size, number of threads in a block" << endl;
      cout << "(2) dimension of the problem" << endl;
      cout << "(3) number of monomials" << endl;
      cout << "(4) number of variables in each monomial" << endl;
      cout << "(5) maximal degree of each variable" << endl;
      cout << "(6) number of repeated runs" << endl;
      cout << "(7) execution mode" << endl;
      cout << "please try again..." << endl; return 1;
   }
   *BS = atoi(argv[1]);     // block size
   *dim = atoi(argv[2]);    // dimension
   *NM = atoi(argv[3]);     // number of monomials
   *NV = atoi(argv[4]);     // number of variables in each monomial
   *deg = atoi(argv[5]);    // maximal degree of each variable
   *r = atoi(argv[6]);      // number of repeated runs
   *mode = atoi(argv[7]);   // execution mode
   return 0;
}

void random_positions ( int dim, int NM, int NV, int *p_int, char *p_char )
// Generates random positions for the NM monomials in NV variables,
// for a system of dimension dim, stored twice in the arrays p_int and p_char.
{
   for(int i=0; i<NM; i++)
   {
      int lp = 0;
      for(int j=0; j<NV; j++)
      {
         int ia = max(rand() % (dim-lp - (NV-j)),1);
         lp = lp + ia;
         p_int[NV*i+j] = (i+j) % dim;
         p_char[NV*i+j] = (char) p_int[NV*i+j];
      }
   }
}

void random_exponents ( int NM, int NV, int d, int *e_int, char *e_char )
// Generates random exponents for the NM monomials in NV variables,
// of degree at most d, stored twice in the arrays e_int and e_char.
{
   for(int i=0; i<NM; i++)
   {
      for(int j=0;j<NV;j++)
      {
         e_int[NV*i+j] = rand() % d;
         e_char[NV*i+j] = (char) e_int[NV*i+j];
      }
   }
}
