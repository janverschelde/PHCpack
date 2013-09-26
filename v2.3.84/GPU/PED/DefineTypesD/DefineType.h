/* This file defines the precision level on host and card.
 * Type T is for the card and is either double, gdd_real, or gqd_real.
 * The corresponding type T1 on the host is double, dd_real, or qd_real.
 * The definition of the precision level is set at compile time with
 * the gcc flag "-D precision=d" for double precision
 *              "-D precision=dd" for double double precision, and
 *              "-D precision=qd" for quad double precision. */

#ifndef __DEFINE_TYPE_H__
#define __DEFINE_TYPE_H__

#include <gqd_type.h>
#include <qd/qd_real.h>

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define p precision
#else
#define p 0
#endif

#if(p == 0)
typedef double T;
typedef double T1;
#elif(p == 1)
typedef gdd_real T;
typedef dd_real T1;
#else
typedef gqd_real T;
typedef qd_real T1;
#endif

#endif
