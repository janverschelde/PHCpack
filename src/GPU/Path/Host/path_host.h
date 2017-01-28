/* path_host.h, created on Dec 6, 2014, by yxc with edits by jv */

#ifndef _PATH_HOST_H_
#define _PATH_HOST_H_

#include "newton_host.h"
#include "predictor_host.h"

template <class ComplexType, class RealType>
bool path_tracker
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter,
   double& timeSec_Predict, double& timeSec_Eval, double& timeSec_MGS,
   int reverse = 0, int verbose = 0 );
/*
 * DESCRIPTION :
 *   Runs the path tracker, tracking one single solution path,
 *   defined by the homotopy in the cpu_inst_hom. */

#include "path_host.tpp"

#endif /* _PATH_HOST_H_ */
