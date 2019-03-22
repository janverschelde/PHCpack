#include <unistd.h>

int corecount ( void )
/* returns the number of available cores */
{
   int number_of_cores = sysconf(_SC_NPROCESSORS_ONLN);

   return number_of_cores;
}
