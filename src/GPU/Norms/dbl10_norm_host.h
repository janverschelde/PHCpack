// The file dbl10_norm_host.h specifies functions to compute
// the 2-norm of a real deca double vector and to normalize a vector.

#ifndef __dbl10_norm_host_h__
#define __dbl10_norm_host_h__

void make_copy
 ( int dim, double *orgrtb, double *orgrix, double *orgrmi, double *orgrrg,
   double *orgrpk, double *orgltb, double *orglix, double *orglmi, 
   double *orglrg, double *orglpk, double *duprtb, double *duprix,
   double *duprmi, double *duprrg, double *duprpk, double *dupltb,
   double *duplix, double *duplmi, double *duplrg, double *duplpk );
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as ten arrays for all ten parts of a deca double.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   orgrtb    highest parts of the original vector of dimension dim;
 *   orgrix    second highest parts of the original vector of dimension dim;
 *   orgrmi    third highest parts of the original vector of dimension dim;
 *   orgrrg    fourth highest parts of the original vector of dimension dim;
 *   orgrpk    fifth highest parts of the original vector of dimension dim;
 *   orgltb    fifth lowest parts of the original vector of dimension dim;
 *   orglix    fourth lowest parts of the original vector of dimension dim;
 *   orglmi    third lowest parts of the original vector of dimension dim;
 *   orglrg    second lowest parts of the original vector of dimension dim;
 *   orglpk    lowest parts of the original vector of dimension dim;
 *   duprtb    space for the highest parts of a vector of dimension dim;
 *   duprix    space for the 2nd highest parts of a vector of dimension dim;
 *   duprmi    space for the 3nd highest parts of a vector of dimension dim;
 *   duprrg    space for the 4th highest parts of a vector of dimension dim;
 *   duprpk    space for the 5th highest parts of a vector of dimension dim;
 *   dupltb    space for the 5th lowest parts of a vector of dimension dim;
 *   duplix    space for the 4th lowest parts of a vector of dimension dim;
 *   duplmi    space for the 3rd lowest parts of a vector of dimension dim;
 *   duplrg    space for the 2nd lowest parts of a vector of dimension dim;
 *   duplpk    space for the lowest parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duprtb    duplicated vector with same contents as orgrtb;
 *   duprix    duplicated vector with same contents as orgrix;
 *   duprmi    duplicated vector with same contents as orgrmi;
 *   duprrg    duplicated vector with same contents as orgrrg;
 *   duprpk    duplicated vector with same contents as orgrpk;
 *   dupltb    duplicated vector with same contents as orgltb;
 *   duplix    duplicated vector with same contents as orglix;
 *   duplmi    duplicated vector with same contents as orglmi;
 *   duplrg    duplicated vector with same contents as orglrg;
 *   duplpk    duplicated vector with same contents as orglpk. */

void CPU_norm
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *normrtb, double *normrix, double *normrmi,
   double *normrrg, double *normrpk, double *normltb, double *normlix,
   double *normlmi, double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v, given by ten arrays,
 *   for all ten parts of a deca double.
 *
 * ON ENTRY :
 *   vrtb      highest parts of a vector of size dim;
 *   vrix      second highest parts of a vector of size dim;
 *   vrmi      third highest parts of a vector of size dim;
 *   vrrg      fourth highest parts of a vector of size dim;
 *   vrpk      fifth highest parts of a vector of size dim;
 *   vltb      fifth lowest parts of a vector of size dim;
 *   vlix      fourth lowest parts of a vector of size dim;
 *   vlmi      third lowest parts of a vector of size dim;
 *   vlrg      second lowest parts of a vector of size dim;
 *   vlpk      lowest parts of a vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   normrtb   highest part of the 2-norm of the vector v;
 *   normrix   second highest part of the 2-norm of the vector v;
 *   normrmi   third highest part of the 2-norm of the vector v;
 *   normrrg   fourth highest lowest part of the 2-norm of the vector v;
 *   normrpk   fifth highest part of the 2-norm of the vector v;
 *   normltb   fifth lowest part of the 2-norm of the vector v;
 *   normlix   fourth lowest part of the 2-norm of the vector v;
 *   normlmi   third lowest part of the 2-norm of the vector v;
 *   normlrg   second lowest lowest part of the 2-norm of the vector v;
 *   normlpk   lowest part of the 2-norm of the vector v. */

void CPU_normalize
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double normrtb, double normrix, double normrmi, double normrrg,
   double normrpk, double normltb, double normlix, double normlmi,
   double normlrg, double normlpk );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vrtb      highest parts of a vector of size dim;
 *   vrix      second highest parts of a vector of size dim;
 *   vrmi      third highest parts of a vector of size dim;
 *   vrrg      fourth highest parts of a vector of size dim;
 *   vrpk      fifth highest parts of a vector of size dim;
 *   vltb      fifth lowest parts of a vector of size dim;
 *   vlix      fourth lowest parts of a vector of size dim;
 *   vlmi      third lowest parts of a vector of size dim;
 *   vlrg      second lowest parts of a vector of size dim;
 *   vlpk      lowest parts of a vector of size dim;
 *   normrtb   highest part of a nonzero number to divide v with;
 *   normrix   second highest part of a nonzero number to divide v with;
 *   normrmi   third highest part of a nonzero number to divide v with;
 *   normrrg   fourth highest part of a nonzero number to divide v with;
 *   normrpk   fifth highest part of a nonzero number to divide v with;
 *   normltb   fifth lowest part of a nonzero number to divide v with.
 *   normlix   fourth lowest part of a nonzero number to divide v with.
 *   normlmi   third lowest part of a nonzero number to divide v with.
 *   normlrg   second lowest part of a nonzero number to divide v with.
 *   normlpk   lowest part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vrtb      highest parts of the normalized vector;
 *   vrix      second highest parts of the normalized vector;
 *   vrmi      third highest parts of the normalized vector;
 *   vrrg      fourth highest parts of the normalized vector;
 *   vrpk      fifth highest parts of the normalized vector;
 *   vltb      fifth lowest parts of the normalized vector;
 *   vlix      fourth lowest parts of the normalized vector;
 *   vlmi      third lowest parts of the normalized vector;
 *   vlrg      second lowest parts of the normalized vector;
 *   vlpk      lowest parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
