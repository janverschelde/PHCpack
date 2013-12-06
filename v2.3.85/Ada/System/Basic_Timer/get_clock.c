/* calls the C routine clock */

#include <time.h>

long get_clock(void)
{
   return clock();
}
