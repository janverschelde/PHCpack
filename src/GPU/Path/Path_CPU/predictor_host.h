/*
 * predictor.h
 *
 *  Created on: Dec 6, 2014
 *      Author: yxc
 */

#ifndef PREDICTOR_HOST_H_
#define PREDICTOR_HOST_H_

#include "DefineType_Host.h"

void predictor_newton(CT** x_array, CT* t_array_orig, int x_t_idx, int np_predictor, int dim);


#endif /* PREDICTOR_HOST_H_ */
