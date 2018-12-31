/* file timer.h contains the prototype for functions to time c code */

#include <time.h>

typedef struct timedata timer;

struct timedata
{
   clock_t user_time;
   time_t  time1, time2;
};

void tstart ( timer *tmr );
/* initializes the timer */

void tstop ( timer *tmr );
/* stop the timer */

void tprint ( timer tmr );
/* formatted print of user CPU time*/
