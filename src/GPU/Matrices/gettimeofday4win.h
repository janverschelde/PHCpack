// The file wingettimeofday.h defines code for the gettimeofday() function 
// to measure wall clock time on Windows.

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

int gettimeofday ( struct timeval * tp, struct timezone * tzp );
/*
 * DESCRIPTION :
 *   Returns the current time in tp. */

#endif
