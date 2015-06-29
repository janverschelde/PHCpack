/*
 * DefineType_Host.h
 *
 *  Created on: Dec 24, 2014
 *      Author: yxc
 */

#ifndef DEFINE_TYPE_HOST_H_
#define DEFINE_TYPE_HOST_H_

#include "../Complex/complexH.h"

typedef double T1;

#define CT complexH<double>

inline double read_number(const char* number){
	return atof(number);
}


#endif /* DEFINE_TYPE_HOST_H_ */
