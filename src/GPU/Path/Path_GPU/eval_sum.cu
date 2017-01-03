#include "eval_sum_seq.cu"
#include "eval_sum_tree.cu"

void eval_sum ( GPUWorkspace& workspace, const GPUInst& inst )
{
   int sum_method = 1;
   if(sum_method == 0)
   {
      eval_sum_seq(workspace, inst);
   }
   else
   {
      eval_sum_tree(workspace, inst);
   }
}
