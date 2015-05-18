/* Solves a hypersurface Pieri problem providing encapsulations
   for the use_c2pieri functionality.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#include <stdio.h>
#include <stdlib.h>

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2pieri ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2pieri ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int initialize_dimensions ( int m, int p, int q );
/*
 * Initializes the dimensions for the Pieri problem:
 *   m : the dimension of the input planes,
 *   p : the dimension of the output planes,
 *   q : the degree of the maps.
 * The dimension of the ambient space is n = m+p.
 * The Pieri homotopies compute all maps of degree q
 * that produce p-planes that meet given m-planes at
 * specified interpolation points.  */

int random_input_planes ( int m, int p, int q );
/*
 * Generates random coefficients of the m*p + q*(m+p)
 * input planes of dimension m and stores these planes
 * into the Pieri homotopy state machine. */

int random_interpolation_points ( int m, int p, int q );
/*
 * Generates random interpolation points to compute
 * all maps of degree that produce p-planes meeting
 * given m-planes at prescribed interpolation points. */

int main ( int argc, char *argv[] )
{
   int m,p,q,fail;

   adainit();

   printf("\nCalling the Pieri homotopies in PHCpack...\n");
   printf("  give m (dimension of input planes) : "); scanf("%d",&m);
   printf("  give p (dimension of output planes) : "); scanf("%d",&p);
   printf("  give q (degree of the maps) : "); scanf("%d",&q);
   fail = initialize_dimensions(m,p,q);
   fail = random_input_planes(m,p,q);

   adafinal();

   return 0;
}

int initialize_dimensions ( int m, int p, int q )
{
   int fail,mpq[3];
   int *b;
   double *c;

   mpq[0] = m; mpq[1] = p; mpq[2] = q;

   fail = _ada_use_c2pieri(1,mpq,b,c);

   return fail;
}

int random_input_planes ( int m, int p, int q )
{
   int n = m*p + q*(m+p);
   int nbcff = 2*(m+p)*m*n;
   int i,fail,mpq[3];
   double c[nbcff];

   mpq[0] = m; mpq[1] = p; mpq[2] = q;     

   for(i=0; i<nbcff; i++)
      c[i] = ((double) rand())/RAND_MAX;

   fail = _ada_use_c2pieri(2,mpq,&n,c);

   return fail;
}

int random_interpolation_points ( int m, int p, int q )
{
   int n = m*p + q*(m+p);
   int nbcff = 2*n;
   int i,fail,mpq[3];
   double c[nbcff];

   mpq[0] = m; mpq[1] = p; mpq[2] = q;
   nbcff = 2*n;

   for(i=0; i<nbcff; i++)
      c[i] = ((double) rand())/RAND_MAX;

   fail = _ada_use_c2pieri(3,mpq,&n,c);
}
