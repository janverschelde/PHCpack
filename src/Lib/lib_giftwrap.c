/* simple test on gift wrapping for convex hulls */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "giftwrappers.h"

void generate_points ( int dim, int m, char *strpts );
/*
 * DESCRIPTION :
 *   Generates m random integer points in dimension dim.
 *   The calling routine must allocate sufficient space
 *   for the array strpts. */

int planar_hull ( int nbc_pts, char *pts );
/*
 * DESCRIPTION :
 *   Given in pts is the list of tuples string representation of
 *   a planar point configuration, with as many characters as nbc_pts.
 *   The giftwrapping method is called and the result is printed. */

int hull_in_3d ( int nbc_pts, char *pts );
/*
 * DESCRIPTION :
 *   Given in pts is the list of tuples string representation of a point
 *   configuration in 3-space, with as many characters as nbc_pts.
 *   The giftwrapping method is called and the result is printed. */

int hull_in_4d ( int nbc_pts, char *pts );
/*
 * DESCRIPTION :
 *   Given in pts is the list of tuples string representation of a point
 *   configuration in 4-space, with as many characters as nbc_pts.
 *   The giftwrapping method is called and the result is printed. */

int main ( void )
{
   int dim,npts;

   adainit();

   printf("Give the ambient dimension : ");
   scanf("%d",&dim);
   printf("Give the number of points : ");
   scanf("%d",&npts);

   {
      char points[dim*10*npts];

      generate_points(dim,npts,points);
      printf("List of tuples : %s\n",points);
      printf("Number of characters : %d\n",strlen(points));
 
      if(dim == 2)
         planar_hull(strlen(points),points);
      else
         if(dim == 3)
            hull_in_3d(strlen(points),points);
         else 
            hull_in_4d(strlen(points),points);
   }

   adafinal();

   return 0;
}

void generate_points ( int dim, int m, char *strpts )
{
   int coords[dim*m],i,j;

   srand(time(NULL));
   for(i=0; i<dim*m; i++)
      coords[i] = -9 + (rand() % 19);  // single digit numbers

   strpts[0] = '[';
   strpts[1] = '\0';
   printf("The coordinates :");
   for(i=0; i<m; i++)
   {
      char strnum[80];
      printf(" %d", coords[i]); 
      sprintf(strnum,"(%d, ", coords[i]);
      strcat(strpts,strnum);
      for(j=1; j<dim-1; j++)
      {
         printf(" %d", coords[j*m+i]); 
         sprintf(strnum,"%d, ", coords[j*m+i]);
         strcat(strpts,strnum);
      }
      printf(" %d", coords[(dim-1)*m+i]);
      if(i<m-1)
         sprintf(strnum,"%d), ", coords[(dim-1)*m+i]);
      else
         sprintf(strnum,"%d)]", coords[(dim-1)*m+i]);
      strcat(strpts,strnum);
   }
   printf("\n");
}

int planar_hull ( int nbc_pts, char *pts )
{
   int fail,nbc_hull;
   char hull[10*nbc_pts];

   fail = convex_hull_2d(nbc_pts,pts,&nbc_hull,hull);

   printf("Vertices and inner normals of the convex hull :\n");
   hull[nbc_hull] = '\0';
   printf("%s\n",hull);
}

int hull_in_3d ( int nbc_pts, char *pts )
{
   int fail,nbr,fcn,ncp;
   char rep[256];

   fail = convex_hull(nbc_pts,pts);
   fail = number_of_facets(3,&nbr);
   printf("Number of facets : %d\n",nbr); 

   do
   {
      printf("Give a facet number (-1 to exit) : ");
      scanf("%d",&fcn);
      if (fcn < 0) break;
      fail = retrieve_facet(3,fcn,&ncp,rep);
      printf("The representation of facet %d :\n %s\n",fcn,rep);
   }
   while (fcn >= 0);

   return fail;
}

int hull_in_4d ( int nbc_pts, char *pts )
{
   int fail,nbr,fcn,ncp;
   char rep[256];

   fail = convex_hull(nbc_pts,pts);
   fail = number_of_facets(4,&nbr);
   printf("Number of facets : %d\n",nbr); 

   do
   {
      printf("Give a facet number (-1 to exit) : ");
      scanf("%d",&fcn);
      if (fcn < 0) break;
      fail = retrieve_facet(4,fcn,&ncp,rep);
      printf("The representation of facet %d :\n %s\n",fcn,rep);
   }
   while (fcn >= 0);

   return fail;
}
