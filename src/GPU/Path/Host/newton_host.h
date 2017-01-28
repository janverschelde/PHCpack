/* newton_host.h, created on: Dec 6, 2014 by yxc with edits by jv */

#ifndef _NEWTON_HOST_H_
#define _NEWTON_HOST_H_

#include "mgs_host.h"
#include "eval_host.h"
#include "parameter.h"

template <class ComplexType, class RealType>
bool CPU_Newton
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse = 0 );

template <class ComplexType, class RealType>
bool CPU_Newton_Refine
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse = 0 );

#include "newton_host.tpp"

#endif /* _NEWTON_HOST_H_ */
