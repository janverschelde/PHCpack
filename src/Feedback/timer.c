
/* file "timer.c" contains definitions of prototypes in "timer.h" */

#include <stdio.h>
#include "timer.h" 
#define MAX_LONG 2147483647


void tstart ( timer *tmr )
{
  tmr->user_time = clock();
  tmr->time1 = time(NULL);
}

void tstop ( timer * tmr )
{
  tmr->user_time = clock() - tmr->user_time;
  tmr->time2 = time(NULL);
}

void tprint ( timer tmr )
{
  int uhour, umin, usec, umsc, remainder,
      thour, tmin, tsec;
  long useconds, tseconds;  

  /* compute the total time */
  tseconds = difftime(tmr.time2, tmr.time1);
  thour = tseconds/3600;
  remainder = tseconds-thour*3600;
  tmin = remainder/60;
  tsec = remainder-tmin*60;
  printf("%s%5d%s%3d%s%3d%s\n",
         "Total elapsed time(resolution of 1 second) :\n",
          thour," hour(s)", tmin," minute(s)",tsec, " second(s).");


  /* compute the user CPU time */
  if(tseconds>(MAX_LONG/CLOCKS_PER_SEC))  
  /* Since the type of clock_t is long integer, it 
     can not hold a number larger than 2,147,483,647 */ 
    {
       printf("User CPU time is not available because of overflow.\n");  
       return; 
    }
  useconds = tmr.user_time/CLOCKS_PER_SEC;
  uhour = useconds/3600;
  remainder = useconds-uhour*3600;
  umin = remainder/60;
  remainder = remainder-umin*60;
  usec = remainder;
  umsc = tmr.user_time*1000/CLOCKS_PER_SEC-useconds*1000;

  printf("%s%5d%s%3d%s%3d%s%4d%s\n",
	 "Elapsed user CPU time :\n", 
          uhour," hour(s)", umin," minute(s)",usec, " second(s)",
	  umsc," millisecond(s).");
} 
