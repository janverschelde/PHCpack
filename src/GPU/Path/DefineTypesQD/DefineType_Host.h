/*
 * DefineType_Host.h
 *
 *  Created on: Dec 24, 2014
 *      Author: yxc
 */

#ifndef __DEFINE_TYPE_QD_HOST_H_
#define __DEFINE_TYPE_QD_HOST_H_

#include "../Complex/complexH.h"

typedef qd_real T1;

#define CT complexH<qd_real>

inline qd_real read_number(const char* number){
	return qd_real(number);
}

#endif /* __DEFINE_TYPE_QD_HOST_H_ */
