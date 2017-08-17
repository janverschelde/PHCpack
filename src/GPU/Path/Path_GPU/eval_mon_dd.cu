#include "eval_mon_single.cu"
#include "eval_mon_seq.cu"
#include "eval_mon_seq_align.cu"
#include "eval_mon_tree.cu"

void eval_mon_tree ( GPUWorkspace& workspace, const GPUInst& inst );

void eval_mon ( GPUWorkspace& workspace, const GPUInst& inst )
{
   if(MON_EVAL_METHOD == 0)
   {
      eval_mon_seq(workspace, inst);
   }
   else if(MON_EVAL_METHOD == 1)
   {
      eval_mon_seq_align(workspace, inst);
   }
   else
   {
      eval_mon_tree(workspace, inst);
   }
}

void eval_mon_tree ( GPUWorkspace& workspace, const GPUInst& inst )
{
   int max_level = 3;

   int* pos_start_tmp = inst.mon_pos_start;
   GT* workspace_coef_tmp = workspace.coef;

   int last_level = min(inst.level, max_level+1);
   for(int i=0; i<last_level; i++)
   {
      if(i==0)
      {
         eval_mon_single_kernel<<<inst.mon_level_grid[0], inst.mon_level0_BS>>>
            (workspace.mon, workspace.x, workspace.coef, inst.mon_pos_start,
             inst.mon_pos, inst.n_mon_level[0]);
      }
      if(i==1)
      {
         eval_mon_tree_2_kernel<<<inst.mon_level_grid[1], inst.mon_global_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp, pos_start_tmp,
             inst.mon_pos, inst.n_mon_level[1]);
      }
      if(i==2)
      {
         eval_mon_tree_4_kernel<<<inst.mon_level_grid[2], inst.mon_level_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp,
             pos_start_tmp, inst.mon_pos, inst.n_mon_level[i]);
      }
      if(i==3)
      {
         eval_mon_tree_kernel<4><<<inst.mon_level_grid[3],inst.mon_level_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp, pos_start_tmp,
             inst.mon_pos, inst.n_mon_level[i]);
      }
      pos_start_tmp += inst.n_mon_level[i];
      workspace_coef_tmp += inst.n_mon_level[i];
   }
   if(inst.level > max_level+1)
   {
      int n_mon_tmp = inst.n_mon_level_rest[max_level];
      if(max_level == 1)
      {
         dim3 mon_level_rest_grid
            = get_grid(n_mon_tmp, inst.mon_global_BS, 1, 1);
         eval_mon_seq_kernel<<<mon_level_rest_grid, inst.mon_global_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp,
             pos_start_tmp, inst.mon_pos, n_mon_tmp);
      }
      else if(max_level == 2)
      {
         dim3 mon_level_rest_grid
            = get_grid(n_mon_tmp, inst.mon_level_BS, 1, 2);
         eval_mon_tree_4n_kernel<<<mon_level_rest_grid, inst.mon_level_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp,
             pos_start_tmp, inst.mon_pos, n_mon_tmp);
      }
      else if(max_level == 3)
      {
         dim3 mon_level_rest_grid
            = get_grid(n_mon_tmp, inst.mon_level_BS, 1, 4);
         eval_mon_tree_n_kernel<4><<<mon_level_rest_grid, inst.mon_level_BS>>>
            (workspace.mon, workspace.x, workspace_coef_tmp,
             pos_start_tmp, inst.mon_pos, n_mon_tmp);
      }
   }
}
