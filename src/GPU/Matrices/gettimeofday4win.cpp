// The file gettimeofday4win.cpp defines the function specified in
// the file gettimeofday4win.h.

#include "gettimeofday4win.h"

int gettimeofday ( struct timeval * tp, struct timezone * tzp )
{
   static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

   SYSTEMTIME  system_time;
   FILETIME    file_time;
   uint64_t    time;

   GetSystemTime( &system_time );
   SystemTimeToFileTime( &system_time, &file_time );
   time =  ((uint64_t)file_time.dwLowDateTime )      ;
   time += ((uint64_t)file_time.dwHighDateTime) << 32;

   tp->tv_sec  = (long) ((time - EPOCH) / 10000000L);
   tp->tv_usec = (long) (system_time.wMilliseconds * 1000);

   return 0;
}
