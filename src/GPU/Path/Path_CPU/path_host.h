/* path_host.h, created on Dec 6, 2014, by yxc with edits by jv */

#ifndef PATH_HOST_H_
#define PATH_HOST_H_

#include "newton_host.h"
#include "predictor_host.h"

bool path_tracker
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Predict, double& timeSec_Eval,
   double& timeSec_MGS, int reverse = 0, int verbose = 0 );
/*
 * DESCRIPTION :
 *   Runs the path tracker, tracking one single solution path, defined
 *   by the homotopy in the cpu_inst_hom. */

#endif /* PATH_HOST_H_ */
