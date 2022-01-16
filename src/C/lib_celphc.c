/* Simple test program in C for the powers of t in a polyhedral homotopy. */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"
#include "syscon.h"
#include "celcon.h"

int compute_mixed_volume ( void );
/*
 * DESCRIPTION :
 *   prompts the user for a polynomial system and then computes
 *   its mixed volume. */

int write_lifted_supports ( int n );
/* 
 * DESCRIPTION :
 *   writes the lifted supports, the points are n-vectors to screen,
 *   returns 0 if no failure occurred, otherwise 1 is returned. */

double inner_product ( int n, double *p, double *q );
/*
 * DESCRIPTION :
 *   Returns the inner product of two n-dimensional arrays
 *   of doubles, stored in p and q respectively. */

int write_inner_products ( int n, int r, int k );
/*
 * DESCRIPTION :
 *    For each point in the lifted support, writes its inner product
 *    with the inner normal to the k-th mixed cell.
 *    The ambient dimension is n and r is the number of different supports. */

int query_cell ( int n, int r, int *cellnb, int tosolve );
/*
 * DESCRIPTION :
 *   Prompts the user for a cell number and lists the cell,
 *   where r is the number of different supports;
 *   returns 0 if no failure occurred, otherwise 1 is returned.
 *   If tosolve = 1, then the start system corresponding to the cell is
 *   solve, otherwise (if tosolve /= 1), there is no solving.
 *   The number of the cell that was queried is in cellnb on return. */

int show_mixture ( int dim, int *r );
/* 
 * DESCRIPTION :
 *   queries the cells container for the type of mixture,
 *   the number "dim" on input should be at least as large as the
 *   number of different supports, if the query results in no failure,
 *  then the mixture is shown and *r contains #different supports. */

void read_and_retrieve ( void );
/*
 * DESCRIPTION :
 *   reads a mixed-cell configuration, initializes the cells container
 *   and tests the retrieval operations. */

int main ( void )
{
   printf("\nComputing the powers of t in a polyhedral homotopy...\n");

   adainit();

   compute_mixed_volume();
   read_and_retrieve();

   adafinal();

   return 0;
}

int compute_mixed_volume( void )
{
   int fail,mv,len,dim,r;

   fail = syscon_read_standard_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_standard_system();
   fail = mixed_volume(&mv);
   printf("\nThe mixed volume : %d\n",mv);
   fail = celcon_number_of_cells(&len);
   printf("\nnumber of mixed cells : %d\n",len);
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);
   fail = show_mixture(dim,&r);
}

int write_lifted_supports ( int n )
{
   int fail,a,b,r,i,j,k;
   int nl[n];
   double *c,pt[n];
 
   fail = celcon_length_of_supports(&r,nl);
   printf("\nlength of the lifted support lists :");
   for(i=0; i<r; i++) printf(" %d",nl[i]); printf("\n");
   printf("the lifted supports:\n");
   for(i=0; i<r; i++)
   {
      for(j=0; j<nl[i]; j++)
      {
        /* printf("retrieving point %d from list %d ...\n",j+1,i+1); */
         fail = celcon_get_lifted_point(n,i+1,j+1,pt);
         for(k=0; k<n-1; k++) printf(" %d",(int) pt[k]);
         printf(" %.15le\n", pt[n-1]);
      }
      printf("\n");
   }
   return fail;
}

double inner_product ( int n, double *p, double *q )
{
   double result = 0.0;
   int i;

   for(i=0; i<n; i++)
      result = result + p[i]*q[i];

   return result;
}

int write_inner_products ( int n, int r, int k )
{
   int fail,i,j,kk;
   int nl[r];
   double normal[n],point[n],prod,minprod;
   char ch;

   fail = celcon_get_inner_normal(n,k,normal);
   if (fail==1)
      printf("an error happened...\n");
   else
   {
      printf("inner normal for cell %d :\n", k);
      for(i=0;i<n;i++) printf("%.15le\n", normal[i]);
      printf("Hit enter to continue.");
      scanf("%c",&ch);
      printf("\ncomputing the inner products with all points ...\n");
      fail = celcon_length_of_supports(&r,nl);
      for(i=0;i<r;i++)
      {
         printf("support %d :\n", i);
         for(j=0;j<nl[i];j++)
         {
            fail = celcon_get_lifted_point(n,i+1,j+1,point);
            for(kk=0; kk<n-1; kk++) printf(" %d",(int) point[kk]);
            prod = inner_product(n,normal,point);
            printf(" : %.15le\n", prod);
            if(j == 0)
               minprod = prod;
            else
               if(prod < minprod) minprod = prod;
         }
         printf("-> support %d has minimal product : %.15le\n", i, minprod);
      }
   }
   return fail;
}

int query_cell ( int n, int r, int *cellnb, int tosolve )
{
   int k,fail,*b,i,j,mv,nl[r],kk,mvcell;
   double *c,normal[n],pt[n];
   char ch;

   printf("Give a cell number : "); scanf("%d", &k);

   if(tosolve == 1)
   {
      fail = celcon_standard_random_coefficient_system();
      fail = celcon_standard_polyhedral_homotopy();
   }

   fail = celcon_get_inner_normal(n,k,normal);
   if (fail==1)
      printf("an error happened...\n");
   else
   {
      *cellnb = k;
      printf("inner normal for cell %d :\n", k);
      for(i=0;i<n;i++) printf("%.15le\n", normal[i]);
      fail = celcon_number_of_points_in_cell(k,r,nl);
      printf("number of points in supports :");
      for(i=0;i<r;i++) printf(" %d",nl[i]); printf("\n");
      scanf("%c",&ch); /* get previous new line symbol */
      printf("Hit enter to continue.");
      scanf("%c",&ch); /* catch new line symbol */
      printf("points in the supports :\n");
      for(i=0;i<r;i++)
      {
         for(j=0;j<nl[i];j++)
         {
            fail = celcon_get_point_in_cell(n,k,i+1,j+1,pt);
            for(kk=0; kk<n-1; kk++) printf(" %d",(int) pt[kk]);
            printf(" %.15le\n", pt[n-1]);
         }
      }
      fail = celcon_mixed_volume(k,&mv);
      printf("mixed volume : %d\n",mv);
      printf("Hit enter to continue.");
      scanf("%c",&ch); /* catch new line symbol */
      if(tosolve == 1)
      {
         fail = celcon_solve_standard_start_system(k,&mvcell);
         printf("The number of solutions : %d\n", mvcell);
         printf("Hit enter to continue.");
         scanf("%c",&ch); /* catch new line symbol */
      }
      {
         int cl[1+r+2*n];
         double inner_normal[n];

         fail = celcon_retrieve_mixed_cell(n,r,k,cl,inner_normal);
         printf("the inner normal for cell %d : \n", k);
         for(i=0; i<n; i++) printf(" %.15le\n", inner_normal[i]);
         printf("total number of points in supports : %d\n",cl[0]);
         printf("number of points in supports : ");
         for(i=1; i<=r; i++) printf(" %d", cl[i]); printf("\n");
         printf("labels of points : ");
         kk=r;
         for(i=1; i<=r; i++)
         {
            for(j=0; j<cl[i]; j++) printf(" %d", cl[++kk]);
            printf(" | ");
         }
         printf("\n");
         printf("Hit enter to continue.");
         scanf("%c",&ch); /* catch new line symbol */
      }
   }
   return fail;
}

int show_mixture ( int dim, int *r )
{
   int fail,*mix,i;
   double *c;

   mix = (int*)calloc(dim,sizeof(int));
   fail = celcon_type_of_mixture(r,mix);

   if(fail == 1)
      printf("An error occurred, type of mixture not available.\n");
   else
   {
      printf("number of different supports : %d\n",*r);
      printf("type of mixture :");
      for(i=0; i<*r; i++) printf(" %d", mix[i]); printf("\n");
   }

   return fail;
}

void read_and_retrieve ( void )
{
   int n,fail,*d,r;
   double *c;
   int len,dim,cell;
   char ans;

   printf("\nDo you wish to see the cells ? (y/n) ");
   ans = getchar();
   if(ans == 'y')
      fail = celcon_write_mixed_cell_configuration();

   fail = celcon_number_of_cells(&len);
   printf("\nnumber of mixed cells : %d\n",len);
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);

   fail = show_mixture(dim,&r);

   fail = write_lifted_supports(dim);
   fail = query_cell(dim,r,&cell,1);
   fail = write_inner_products(dim,r,cell);
}
