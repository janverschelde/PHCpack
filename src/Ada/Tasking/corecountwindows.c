#include <windows.h>

int corecount ( void )
/* returns the number of available cores on Windows */
{
   SYSTEM_INFO sysinfo;
   GetSystemInfo(&sysinfo);
   int number_of_cores = sysinfo.dwNumberOfProcessors;

   return number_of_cores;
}
