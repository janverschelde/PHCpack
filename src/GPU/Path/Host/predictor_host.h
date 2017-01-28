/* predictor.h, created on Dec 6, 2014 by yxc with edits by jv */

#ifndef _PREDICTOR_HOST_H_
#define _PREDICTOR_HOST_H_

template <class ComplexType>
void predictor_newton
 ( ComplexType** x_array, ComplexType* t_array_orig,
   int x_t_idx, int np_predictor, int dim );
/*
 * DESCRIPTION :
 *   Applies Newton extrapolation to predict the next point
 *   on a solution path. */

#include "predictor_host.tpp"

#endif /* _PREDICTOR_HOST_H_ */
