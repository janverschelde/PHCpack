// The file dbl_onenorms_host.h specifies functions to compute the 1-norm
// of vectors in double precision.

#ifndef __dbl_onenorms_host_h__
#define __dbl_onenorms_host_h__

void CPU_dbl_onenorm ( int dim, double *v, double *nrm );
/*
 * DESCRIPTION :
 *   Returns in nrm the 1-norm of the vector v of dimension dim,
 *   which is the largest absolute value of the elements in v. */

void CPU_cmplx_onenorm ( int dim, double *vre, double *vim, double *nrm );
/*
 * DESCRIPTION :
 *   Returns in nrm the 1-norm of the vector v of dimension dim,
 *   which is the largest sum of the absolute value of real and
 *   imaginary parts of the elements in v.
 *
 * ON ENTRY :
 *   dim     dimension of the vector v;
 *   vre     real parts of the elements in v;
 *   vim     imaginary parts of the elements in v. */

#endif
