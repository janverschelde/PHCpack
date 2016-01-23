/* simple test on managing numerically computed tropisms */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "numbtrop.h"

void generate_winding_numbers ( int nbt, int *wind );
/*
 * DESCRIPTION :
 *   Generates as many random single digit integers as the value of nbt.
 *   The random numbers are stored in wind. */

void random_doubles ( int n, double *d );
/*
 * DESCRIPTION :
 *   Generates n random doubles and stores the numbers in d. */

void generate_directions ( int nbt, int dim, double *dirs );
/*
 * DESCRIPTION :
 *   Generates nbt*dim random doubles and returns those in dirs. */

void generate_errors ( int nbt, double *err );
/*
 * DESCRIPTION :
 *   Generates as many random doubles as the value of nbt,
 *   and stores the values in err. */

void standard_generate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Generates random numbers for nbt tropisms of length dim,
 *   with coordinates stored in dir, winding numbers in wind,
 *   and errors in err, in standard double precision.
 *   Data must have been allocated for wind, dir, and err. */

void standard_operate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Tests the operations for nbt tropisms of length dim,
 *   with coordinates in dir, winding numbers in wind, and
 *   errors in err, in standard double precision. */

void standard_test ( int nbt, int dim );
/*
 * DESCRIPTION :
 *   Tests for nbt numerically computed tropisms of dimension dim,
 *   computed in standard double precision. */

int main ( void )
{
   int dim,nbt;

   adainit();

   printf("Give the dimension : "); scanf("%d",&dim);
   printf("Give the number of directions : "); scanf("%d",&nbt);

   srand(time(NULL));
   standard_test(nbt,dim);

   adafinal();

   return 0;
}

void generate_winding_numbers ( int nbt, int *wind )
{
   int i;

   for(i=0; i<nbt; i++)
      wind[i] = -9 + (rand() % 19);  // single digit numbers
}

void random_doubles ( int n, double *d )
{
   int i;
   for(i=0; i<n; i++)
      d[i] = ((double) rand())/RAND_MAX;
}

void generate_directions ( int nbt, int dim, double *dirs )
{
   random_doubles(nbt*dim,dirs);
}

void generate_errors ( int nbt, double *err )
{
   int i;
   random_doubles(nbt,err);
   for(i=0; i<nbt; i++) err[i] = err[i]/1000.0;
}

void standard_generate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int i,j,k;

   generate_winding_numbers(nbt,wind);
   printf("The winding numbers :");
   for(k=0; k<nbt; k++) printf(" %d", wind[k]);
   printf("\n");
   generate_directions(nbt,dim,dir);
   printf("The directions :\n");
   k = 0;
   for(i=0; i<nbt; i++)
   {
      printf("  direction %d :",i+1);
      for(j=0; j<dim; j++) printf(" %.3e",dir[k++]);
      printf("\n");
   }
   generate_errors(nbt,err);
   printf("The errors :\n");
   for(k=0; k<nbt; k++) printf(" %.3e",err[k]);
   printf("\n");
}

void standard_operate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int fail,size,k,idx,wnd;

   fail = standard_initialize(nbt,dim,wind,dir,err);
   fail = standard_size(&size);
   printf("The number of tropisms : %d\n",size);

   do
   {
      double diridx[dim];
      double erridx;

      printf("Give an index (<= 0 to exit): "); scanf("%d",&idx);
      if(idx <= 0) break;
      fail = standard_retrieve_tropism(dim,idx,&wnd,diridx,&erridx);
      printf("The tropism %d has winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e\n",erridx);

      generate_winding_numbers(1,&wnd);
      random_doubles(dim,diridx);
      random_doubles(1,&erridx);

      printf("Changing tropism %d to winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e\n",erridx);
      fail = store_standard_tropism(dim,idx,wnd,diridx,erridx);
   }
   while(idx > 0);

   fail = standard_clear();
}

void standard_test ( int nbt, int dim )
{
   int wind[nbt];
   double dir[nbt*dim];
   double err[nbt];

   standard_generate(nbt,dim,wind,dir,err);
   standard_operate(nbt,dim,wind,dir,err);
}
