/* This file defines the precision level on host and card.
 * Type T is for the card and is either double, gdd_real, or gqd_real.
 * The corresponding type T1 on the host is double, dd_real, or qd_real.
 * The definition of the precision level is set at compile time with
 * the gcc flag "-D precision=d" for double precision
 *              "-D precision=dd" for double double precision, and
 *              "-D precision=qd" for quad double precision. */

#ifndef __DEFINE_TYPE_DD_H__
#define __DEFINE_TYPE_DD_H__

#include <gqd_type.h>
#include <qd/qd_real.h>

typedef gqd_real T;
typedef qd_real T1;

#endif
