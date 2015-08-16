/* predictor.h, created on Dec 6, 2014 by yxc with edits by jv */

#ifndef PREDICTOR_HOST_H_
#define PREDICTOR_HOST_H_

#include "DefineType_Host.h"

void predictor_newton
 ( CT** x_array, CT* t_array_orig, int x_t_idx, int np_predictor, int dim );
/*
 * DESCRIPTION :
 *   Applies Newton extrapolation to predict the next point
 *   on a solution path. */

#endif /* PREDICTOR_HOST_H_ */
