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

void standard_test ( int nbt, int dim )
{
   int i,j,k,fail,size;
   int wind[nbt];
   double dir[nbt*dim];
   double err[nbt];

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

   fail = standard_initialize(nbt,dim,wind,dir,err);
   fail = standard_size(&size);
   printf("The number of tropisms : %d\n",size);
   {
      int idx,wnd;
      double diridx[dim];
      double erridx;
      printf("Give an index : "); scanf("%d",&idx);
      fail = standard_retrieve_tropism(dim,idx,&wnd,diridx,&erridx);
      printf("The tropism %d has winding number %d.\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e\n",erridx);
   }
   fail = standard_clear();
}
