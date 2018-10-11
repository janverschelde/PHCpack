// The file predictor_host.h defines the Newton extrapolating predictor
// for the  path tracker to run on the host.

#ifndef __PREDICTOR_HOST_H__
#define __PREDICTOR_HOST_H__

template <class ComplexType>
void predictor_divdif
 ( ComplexType** x_array, ComplexType* t_array_orig,
   int x_t_idx, int np_predictor, int dim,
   ComplexType* div_diff, ComplexType* t_array, ComplexType* t_diff );
/*
 * Applies Newton extrapolation as a predictor,
 * with workspace arrays in div_diff, t_array, and t_diff.
 */

template <class ComplexType>
void predictor_newton
 ( ComplexType** x_array, ComplexType* t_array_orig,
   int x_t_idx, int np_predictor, int dim );
/*
 * DESCRIPTION :
 *   Applies Newton extrapolation to predict the next point on a path.
 *   This function is not thread safe as it allocates and deallocates
 *   memory to hold intermediate results. */

#include "predictor_host.tpp"

#endif /* __PREDICTOR_HOST_H__ */
