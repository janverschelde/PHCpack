/* simple test on gift wrapping for convex hulls */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "giftwrappers.h"

void test2d ( int m );
/*
 * DESCRIPTION :
 *   Generates m random integer points in the plane
 *   and calls the gift wrapping method. */

int call_giftwrapper ( int nbc_pts, char *pts );
/*
 * DESCRIPTION :
 *   Given in pts is the list of tuples string representation of
 *   a planar point configuration, with as many characters as nbc_pts.
 *   The giftwrapping method is called and the result is printed. */

int main(void)
{
   int npts;

   printf("Give the number of points : ");
   scanf("%d",&npts);

   adainit();

   test2d(npts);

   adafinal();

   return 0;
}

void test2d ( int m )
{
   int coords[2*m],i;
   char points[2*10*m];

   srand(time(NULL));
   for(i=0; i<2*m; i++)
      coords[i] = -9 + (rand() % 19);  // single digit numbers

   points[0] = '[';
   points[1] = '\0';
   printf("The coordinates :");
   for(i=0; i<m; i++)
   {
      char strnum[80];
      printf(" %d", coords[i]); 
      sprintf(strnum,"(%d, ", coords[i]);
      strcat(points,strnum);
      printf(" %d", coords[m+i]); 
      if(i<m-1)
         sprintf(strnum,"%d), ", coords[m+i]);
      else
         sprintf(strnum,"%d)]", coords[m+i]);
      strcat(points,strnum);
   }
   printf("\n");
   printf("List of tuples : %s\n", points);
   printf("Number of characters : %d\n", strlen(points));
 
   call_giftwrapper(strlen(points),points);
}

int call_giftwrapper ( int nbc_pts, char *pts )
{
   int fail,nbc_hull;
   char hull[10*nbc_pts];

   fail = convex_hull_2d(nbc_pts,pts,&nbc_hull,hull);

   printf("Vertices and inner normals of the convex hull :\n");
   hull[nbc_hull] = '\0';
   printf("%s\n",hull);
}
