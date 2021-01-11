// The file wingettimeofday.h defines code for the gettimeofday() function 
// to measure wall clockk time on Windows.

#ifndef __wingettimeofday_h__
#define __wingettimeofday_h__

#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 

// MSVC defines timeval in winsock2.h

/* 
typedef struct timeval
{
   long tv_sec;
   long tv_usec;
}
timeval;
*/

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

#endif
