/*
    This file is a component of DEMiCs
    Copyright (C) 2007 Tomohiko Mizutani, Masakazu Kojima and Akiko Takeda

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef __GROBAL_H
#define __GROBAL_H

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>
// #include <sys/times.h>
#include <unistd.h>
#include <cassert>

#ifndef CLK_TCK
#define  CLK_TCK  sysconf(_SC_CLK_TCK)
#endif

using namespace std;

#define ITER_BLAND 25
#define ITER 1000
#define REINV 0

#define EXCESS 5
#define STRLENGTH 128

#define FALSE 0
#define TRUE 1

#define LENGTHOFTOKEN 128
#define PLUSZERO 1.0E-8
#define MINUSZERO -1.0E-8
#define BIGDOUBLE 1.0E+16
#define SMALLDOUBLE -1.0E+16
#define BIGINT 1000000000

#define FREE 2
#define NONNEGATIVE 3
#define POSITIVE 25
#define NEGATIVE 26

#define POSTHETA 6
#define NEGTHETA 7

#define OPT 4
#define UNBOUNDED  8
#define FEASIBLE 10
#define INFEASIBLE 11

#define CONTINUE 9
#define STOP 14
#define PIVOT_IN 28

//#define ART 15
//#define NOART 16
//#define OUTART 20

#define ICHECK 20
#define MCHECK 21

#define TRIANGLE 28
#define SQUARE 29

#define SLIDE 16
#define DOWN 17

#define NODE 22
#define FN 23
#define FNN 24

#define ON 30
#define OFF 31

#define UNB_TAR 32
#define UNB_COR 33

#define ERROR_ITER 27

///// DEBUG /////////////////////////////
#define DBG_INFO 0

#define DBG_NODE 0
#define DBG_FEA 0
#define DBG_SUC 0
#define DBG_TMP 0
#define DEBUG 0

// For relation table
#define DBG_REL_TABLE 0

// For initial one_point_test
#define DBG_INI_CUR_INFO 0

// For one_point_test
#define DBG_PRE_INFO 0
#define DBG_CUR_INFO 0

// For suc_one_point_test
#define DBG_S_PRE_INFO 0
#define DBG_S_CUR_INFO 0

// For solLP
#define DBG_TSOLLP 0
#define DBG_SOLLP 0
#define DBG_ISOLLP 0
#define DBG_MSOLLP 0

// For chooseSup
#define DBG_CHOOSESUP 0
#define DBG_FINDUNB 0

// For modify p_sol
#define DBG_MODIFY 0

#endif
