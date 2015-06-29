
/* This file defines the precision level on host and card.
 * Type T is for the card and is either double, gdd_real, or gqd_real.
 * The corresponding type T1 on the host is double, dd_real, or qd_real.
 * The definition of the precision level is set at compile time with
 * the gcc flag "-D precision=d" for double precision
 *              "-D precision=dd" for double double precision, and
 *              "-D precision=qd" for quad double precision. */

#ifndef __DEFINE_TYPE_DD_HOST_H__
#define __DEFINE_TYPE_DD_HOST_H__

#include "../Complex/complexH.h"

typedef dd_real T1;

#define CT complexH<dd_real>

inline dd_real read_number(const char* number){
	return dd_real(number);
}

#endif
