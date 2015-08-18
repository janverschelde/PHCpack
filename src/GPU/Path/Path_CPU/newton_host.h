/* newton_host.h, created on: Dec 6, 2014 by yxc with edits by jv */

#ifndef NEWTON_HOST_H_
#define NEWTON_HOST_H_

#include "mgs_host.h"
#include "DefineType_Host.h"
#include "eval_host.h"
#include "parameter.h"

bool CPU_Newton
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse = 0 );


bool CPU_Newton_Refine
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse = 0 );

#endif /* NEWTON_HOST_H_ */
