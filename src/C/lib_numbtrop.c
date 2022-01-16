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

void dobldobl_generate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Generates random numbers for nbt tropisms of length dim,
 *   with coordinates stored in dir, winding numbers in wind,
 *   and errors in err, in double double precision.
 *   Data must have been allocated for wind, dir, and err.
 *   Note that dir occupies 2*nbt*dim doubles and err 2*nbt doubles. */

void quaddobl_generate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Generates random numbers for nbt tropisms of length dim,
 *   with coordinates stored in dir, winding numbers in wind,
 *   and errors in err, in quad double precision.
 *   Data must have been allocated for wind, dir, and err.
 *   Note that dir occupies 4*nbt*dim doubles and err 4*nbt doubles. */

void standard_operate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Tests the operations for nbt tropisms of length dim,
 *   with coordinates in dir, winding numbers in wind, and
 *   errors in err, in standard double precision. */

void dobldobl_operate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Tests the operations for nbt tropisms of length dim,
 *   with coordinates in dir, winding numbers in wind, and
 *   errors in err, in double double precision.
 *   Note that dir has 2*nbt*dim doubles and that err
 *   holds space for 2*nbt doubles. */

void quaddobl_operate
 ( int nbt, int dim, int *wind, double *dir, double *err );
/*
 * DESCRIPTION :
 *   Tests the operations for nbt tropisms of length dim,
 *   with coordinates in dir, winding numbers in wind, and
 *   errors in err, in quad double precision.
 *   Note that dir has 4*nbt*dim doubles and that err
 *   holds space for 4*nbt doubles. */

void standard_test ( int nbt, int dim );
/*
 * DESCRIPTION :
 *   Tests for nbt numerically computed tropisms of dimension dim,
 *   stored in standard double precision. */

void dobldobl_test ( int nbt, int dim );
/*
 * DESCRIPTION :
 *   Tests for nbt numerically computed tropisms of dimension dim,
 *   stored in double double precision. */

void quaddobl_test ( int nbt, int dim );
/*
 * DESCRIPTION :
 *   Tests for nbt numerically computed tropisms of dimension dim,
 *   stored in quad double precision. */

int main ( void )
{
   int dim,nbt,prc;

   adainit();

   printf("Give the dimension : "); scanf("%d",&dim);
   printf("Give the number of directions : "); scanf("%d",&nbt);
   printf("give the precision: 0, 1, or 2 for d, dd, or qd : ");
   scanf("%d",&prc); 

   srand(time(NULL));
   if(prc == 0)
      standard_test(nbt,dim);
   else if(prc == 1)
      dobldobl_test(nbt,dim);
   else if(prc == 2)
      quaddobl_test(nbt,dim);
   else
      printf("Wrong answer for precision.\n");

   adafinal();

   return 0;
}

void generate_winding_numbers ( int nbt, int *wind )
{
   int i;

   for(i=0; i<nbt; i++)
      wind[i] = -99 + (rand() % 199);  // double digit numbers
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

void dobldobl_generate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int i,j,k;

   generate_winding_numbers(nbt,wind);
   printf("The winding numbers :");
   for(k=0; k<nbt; k++) printf(" %d", wind[k]);
   printf("\n");
   generate_directions(nbt,2*dim,dir);
   printf("The directions :\n");
   k = 0;
   for(i=0; i<nbt; i++)
   {
      printf("  direction %d :",i+1);
      for(j=0; j<2*dim; j++) printf(" %.3e",dir[k++]);
      printf("\n");
   }
   generate_errors(2*nbt,err);
   printf("The errors :\n");
   for(k=0; k<2*nbt; k++) printf(" %.3e",err[k]);
   printf("\n");
}

void quaddobl_generate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int i,j,k;

   generate_winding_numbers(nbt,wind);
   printf("The winding numbers :");
   for(k=0; k<nbt; k++) printf(" %d", wind[k]);
   printf("\n");
   generate_directions(nbt,4*dim,dir);
   printf("The directions :\n");
   k = 0;
   for(i=0; i<nbt; i++)
   {
      printf("  direction %d :",i+1);
      for(j=0; j<4*dim; j++) printf(" %.3e",dir[k++]);
      printf("\n");
   }
   generate_errors(4*nbt,err);
   printf("The errors :\n");
   for(k=0; k<4*nbt; k++) printf(" %.3e",err[k]);
   printf("\n");
}

void standard_operate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int fail,size,k,idx,wnd,retdim;

   fail = numbtrop_standard_initialize(nbt,dim,wind,dir,err);
   fail = numbtrop_standard_size(&size);
   printf("The number of tropisms : %d\n",size);
   fail = numbtrop_standard_dimension(&retdim);
   printf("The retrieved dimension : %d\n",retdim);

   do
   {
      double diridx[dim];
      double erridx;

      printf("Give an index (<= 0 to exit): "); scanf("%d",&idx);
      if(idx <= 0) break;
      fail = numbtrop_standard_retrieve_tropism(dim,idx,&wnd,diridx,&erridx);
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
      fail = numbtrop_store_standard_tropism(dim,idx,wnd,diridx,erridx);
   }
   while(idx > 0);

   fail = numbtrop_standard_clear();
}

void dobldobl_operate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int fail,size,k,idx,wnd,retdim;

   fail = numbtrop_dobldobl_initialize(nbt,dim,wind,dir,err);
   fail = numbtrop_dobldobl_size(&size);
   printf("The number of tropisms : %d\n",size);
   fail = numbtrop_dobldobl_dimension(&retdim);
   printf("The retrieved dimension : %d\n",retdim);

   do
   {
      double diridx[2*dim];
      double erridx[2];

      printf("Give an index (<= 0 to exit): "); scanf("%d",&idx);
      if(idx <= 0) break;
      fail = numbtrop_dobldobl_retrieve_tropism(dim,idx,&wnd,diridx,erridx);
      printf("The tropism %d has winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<2*dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e %.3e\n",erridx[0],erridx[1]);

      generate_winding_numbers(1,&wnd);
      random_doubles(2*dim,diridx);
      random_doubles(2,erridx);

      printf("Changing tropism %d to winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<2*dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e %.3e\n",erridx[0],erridx[1]);
      fail = numbtrop_store_dobldobl_tropism(dim,idx,wnd,diridx,erridx);
   }
   while(idx > 0);

   fail = numbtrop_dobldobl_clear();
}

void quaddobl_operate
 ( int nbt, int dim, int *wind, double *dir, double *err )
{
   int fail,size,k,idx,wnd,retdim;

   fail = numbtrop_quaddobl_initialize(nbt,dim,wind,dir,err);
   fail = numbtrop_quaddobl_size(&size);
   printf("The number of tropisms : %d\n",size);
   fail = numbtrop_quaddobl_dimension(&retdim);
   printf("The retrieved dimension : %d\n",retdim);

   do
   {
      double diridx[4*dim];
      double erridx[4];

      printf("Give an index (<= 0 to exit): "); scanf("%d",&idx);
      if(idx <= 0) break;
      fail = numbtrop_quaddobl_retrieve_tropism(dim,idx,&wnd,diridx,erridx);
      printf("The tropism %d has winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<4*dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e %.3e %.3e %.3e\n",
             erridx[0],erridx[1],erridx[2],erridx[3]);

      generate_winding_numbers(1,&wnd);
      random_doubles(4*dim,diridx);
      random_doubles(4,erridx);

      printf("Changing tropism %d to winding number %d,\n",idx,wnd);
      printf("with coordinates :");
      for(k=0; k<4*dim; k++) printf(" %.3e",diridx[k]);
      printf("\nand error : %.3e %.3e %.3e %.3e\n",
             erridx[0],erridx[1],erridx[2],erridx[3]);
      fail = numbtrop_store_quaddobl_tropism(dim,idx,wnd,diridx,erridx);
   }
   while(idx > 0);

   fail = numbtrop_quaddobl_clear();
}

void standard_test ( int nbt, int dim )
{
   int wind[nbt];
   double dir[nbt*dim];
   double err[nbt];

   standard_generate(nbt,dim,wind,dir,err);
   standard_operate(nbt,dim,wind,dir,err);
}

void dobldobl_test ( int nbt, int dim )
{
   int wind[nbt];
   double dir[2*nbt*dim];
   double err[2*nbt];

   dobldobl_generate(nbt,dim,wind,dir,err);
   dobldobl_operate(nbt,dim,wind,dir,err);
}

void quaddobl_test ( int nbt, int dim )
{
   int wind[nbt];
   double dir[4*nbt*dim];
   double err[4*nbt];

   quaddobl_generate(nbt,dim,wind,dir,err);
   quaddobl_operate(nbt,dim,wind,dir,err);
}
