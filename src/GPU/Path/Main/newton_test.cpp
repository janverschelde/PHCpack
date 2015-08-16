/* newton_test.cpp, created on Feb 8, 2015 by yxc with edits by jv */

#include "newton_test.h"

T1 newton_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT t )
{
   std::cout << "--------- Newton Test ----------" << std::endl;

   double timeSec_Eval = 0;
   double timeSec_MGS = 0;

   struct timeval start, end;
   long seconds, useconds;
   double timeMS_cpu;
   gettimeofday(&start, NULL);

   bool success_cpu = CPU_Newton(workspace_cpu,cpu_inst_hom,path_parameter,
                                 timeSec_Eval,timeSec_MGS);
   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   timeMS_cpu = seconds*1000 + useconds/1000.0;

   CT* x_gpu;
   bool success_gpu = GPU_Newton(cpu_inst_hom, path_parameter,sol0,t,x_gpu);
   cout << "Path CPU Newton    Time: " << timeMS_cpu   << endl;
   cout << "Path CPU Eval      Time: " << timeSec_Eval << endl;
   cout << "Path CPU MGS       Time: " << timeSec_MGS  << endl;

   std::cout << "success_cpu = " << success_cpu << std::endl;
   std::cout << "success_gpu = " << success_gpu << std::endl;

   std::cout << "--------- Newton Error Check CPU VS GPU----------"
             << std::endl;
   T1 err = err_check_workspace(workspace_cpu.x, x_gpu, cpu_inst_hom.dim);
   std::cout << "err = " << err << std::endl;
   std::cout << " x_cpu[0] = " << workspace_cpu.x[0];
   std::cout << " x_gpu[0] = " << x_gpu[0];
   delete[] x_gpu;
   x_gpu = NULL;

   return err;
}
