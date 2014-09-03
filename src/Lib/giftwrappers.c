/* the file giftwrappers.c contains the definitions of the functions 
 * with prototypes documented in giftwrappers.h */

#include <stdlib.h>
#include <stdio.h>

extern void adainit ( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal ( void );

int convex_hull_2d ( int nc_pts, char *pts, int *nc_hull, char *hull )
{
   int fail,i;
   int points[10*nc_pts];
   double *c;

   for(i=0; i<nc_pts; i++) points[i] = (int) pts[i];

   fail = _ada_use_c2phc(580,&nc_pts,points,c);

   *nc_hull = nc_pts;

   for(i=0; i<nc_pts; i++) hull[i] = (char) points[i];
  
   return fail;
}
