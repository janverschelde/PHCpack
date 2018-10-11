/* newton_host.h, created on: Dec 6, 2014 by yxc with edits by jv */

#ifndef __NEWTON_HOST_H__
#define __NEWTON_HOST_H__

#include "mgs_host.h"
#include "eval_host.h"
#include "parameter.h"

double sum_norm ( double real, double imag );
/*
 * Returns the sum of the absolute values of real and imag,
 * which are considered the real and imaginary part of a complex number.
 */

double max_norm ( complexH<double>* sol, int dim );
/*
 * Returns the largest sum_norm of the first dim numbers in sol.
 */

double max_norm ( complexH<dd_real>* sol, int dim );
/*
 * Returns the largest sum_norm of the first dim numbers in sol.
 */

double max_norm ( complexH<qd_real>* sol, int dim );
/*
 * Returns the largest sum_norm of the first dim numbers in sol.
 */

template <class ComplexType, class RealType>
bool CPU_Newton
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse=0, bool verbose=false );
/*
 * Runs Newton's method on the system defined by cpu_inst_hom,
 * and on the solutions stored in workspace_cpu,
 * with the tolerances defined in path_parameter for along a path.
 */

template <class ComplexType, class RealType>
bool CPU_Newton_Refine
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse=0, bool verbose=false );
/*
 * Runs Newton's method on the system defined by cpu_inst_hom,
 * and on the solutions stored in workspace_cpu,
 * with the tolerances defined in path_parameter at the end of a path.
 */

#include "newton_host.tpp"

#endif /* __NEWTON_HOST_H__ */
