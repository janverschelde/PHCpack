// The file eval_host.h contains the #include statements for the classes
// to evaluate and differentiate polynomial homotopies.

#ifndef __CPU_INSTRUCTION_EVAL_H__
#define __CPU_INSTRUCTION_EVAL_H__

#include <sys/time.h>
#include <unistd.h>
#include "monset.h"
#include "parameter.h"
#include "workspace_host.h"
#include "path_data.h"
#include "polysys.h"

#define warp_size 32

#include "cpuinsthomcoef.h"
#include "cpuinsthommon.h"
#include "cpuinsthommonblock.h"
#include "cpuinsthomsumblock.h"
#include "cpuinsthomsum.h"
#include "cpuinsthomeq.h"
#include "cpuinsthom.h"

#endif /* __CPU_INSTRUCTION_EVAL_H__ */
