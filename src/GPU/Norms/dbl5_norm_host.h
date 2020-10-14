// The file dbl5_norm_host.h specifies functions to compute
// the 2-norm of a real penta double vector and to normalize a vector.

#ifndef __dbl5_norm_host_h__
#define __dbl5_norm_host_h__

void make_copy
 ( int dim,
   double *orgtb, double *orgix, double *orgmi, double *orgrg, double *orgpk,
   double *duptb, double *dupix, double *dupmi, double *duprg, double *duppk );
/*
 * DESCRIPTION :
 *   Makes a copy of the vector in org to the vector dup,
 *   given as five arrays for all five parts of a penta double.
 *
 * ON ENTRY :
 *   dim       dimension of the vectors org and dup;
 *   orgtb     highest parts of the original vector of dimension dim;
 *   orgix     second highest parts of the original vector of dimension dim;
 *   orgmi     middle parts of the original vector of dimension dim;
 *   orgrg     second lowest parts of the original vector of dimension dim;
 *   orgpk     lowest parts of the original vector of dimension dim;
 *   duptb     space for the highest parts of a vector of dimension dim;
 *   dupix     space for the 2nd highest parts of a vector of dimension dim;
 *   dupmi     space for the middle parts of a vector of dimension dim;
 *   duprg     space for the 2nd lowest parts of a vector of dimension dim;
 *   duppk     space for the lowest parts of a vector of dimension dim.
 *
 * ON RETURN :
 *   duptb     duplicated vector with same contents as orgtb;
 *   dupix     duplicated vector with same contents as orgix;
 *   dupmi     duplicated vector with same contents as orgmi;
 *   duprg     duplicated vector with same contents as orgrg;
 *   duppk     duplicated vector with same contents as orgpk. */

void CPU_norm
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk,
   int dim, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk );
/*
 * DESCRIPTION :
 *   Computes the 2-norm of the vector in v, given by five arrays,
 *   for all five parts of a penta double.
 *
 * ON ENTRY :
 *   vtb       highest parts of a vector of size dim;
 *   vix       second highest parts of a vector of size dim;
 *   vmi       middle parts of a vector of size dim;
 *   vrg       second lowest parts of a vector of size dim;
 *   vpk       lowest parts of a vector of size dim;
 *   dim       dimension of the vector v.
 *
 * ON RETURN :
 *   normtb    highest part of the 2-norm of the vector v;
 *   normix    second highest part of the 2-norm of the vector v;
 *   normmi    middle part of the 2-norm of the vector v;
 *   normrg    second lowest part of the 2-norm of the vector v;
 *   normpk    lowest part of the 2-norm of the vector v. */

void CPU_normalize
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double normtb, double normix, double normmi, double normrg, double normpk );
/*
 * DESCRIPTION :
 *   Normalizes the vector in v.
 *
 * ON ENTRY :
 *   vtb       highest parts of a real vector, of size dim;
 *   vix       second highest parts of a real vector, of size dim;
 *   vmi       middle parts of a real vector, of size dim;
 *   vrg       second lowest parts of a real vector, of size dim;
 *   vpk       lowest parts of a real vector, of size dim;
 *   dim       dimension of the vector v;
 *   normtb    highest part of a nonzero number to divide v with;
 *   normix    second highest part of a nonzero number to divide v with;
 *   normmi    middle part of a nonzero number to divide v with;
 *   normrg    second lowest part of a nonzero number to divide v with.
 *   normpk    lowest part of a nonzero number to divide v with.
 *
 * ON RETURN :
 *   vtb       highest parts of the normalized vector;
 *   vix       second highest parts of the normalized vector;
 *   vmi       middle parts of the normalized vector;
 *   vrg       second lowest parts of the normalized vector;
 *   vpk       lowest parts of the normalized vector,
 *             if on entry, norm equals the 2-norm of v, 
 *             then the 2-norm of v on return will be one. */

#endif
