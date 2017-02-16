#ifndef GPU_DATA_H_
#define GPU_DATA_H_

#include "DefineType.h"
#include "gqd_qd_util.h"
#include "eval_host.h"

int get_NB ( int n_job, int BS, int n_thread_per_job = 1 );

dim3 get_grid ( int NB, int n_path=1 );

dim3 get_grid ( int n_job, int BS, int n_path, int n_thread_per_job = 1 );

class GPUWorkspace
{
   public:

      // Array for predictor usage
      int dim;
      int n_eq;
      int n_array;
      int n_predictor;
      int mon_pos_size;
      int n_constant;
      int n_coef;

      int n_row;
      int n_col;
      int n_matrix;
      int n_matrix_R;
      int n_matrix_P;

      int n_path;
      int n_path_continuous;

      GT* x;
      GT* x_last;
      GT* t;
      GT* t_last;
      int x_t_idx;

      CT alpha;
      GT* alpha_gpu;

      // int arrays
      int* int_arrays;
      int* x_t_idx_mult;
      int* n_point_mult;
      int* path_idx;
      int* path_success;
      int* n_success;
      int* end_range;
      int* newton_success;

      double* double_arrays;
      double* max_delta_x_gpu;
      double* r_max_delta_x_gpu;
      double* max_f_val_gpu;
      double* max_f_val_last_gpu;
      double* r_max_f_val_gpu;
      double* max_x_gpu;

      size_t size_GT_arrays;
      int n_GT_arrays;
      GT* GT_arrays;
      GT* mon;
      GT* coef;
      GT* sum;
      GT* matrix;
      GT* f_val;
      GT* V;
      GT* R;
      GT* P;
      GT* sol;
      GT* matrix_horizontal_mult;
      GT* x_array;
      GT* t_array;
      GT* one_minor_t;

      GT* x_mult;
      GT* t_mult;
      GT* t_last_mult;
      GT* delta_t_mult;
      GT* newton_t_mult;

      int* newton_success_host;
      int* path_idx_host;
      int* path_success_host;

      double* max_delta_x_host;
      double* r_max_delta_x_host;
      double* max_f_val_host;
      double* r_max_f_val_host;
      double* max_x_host;

      GT* workspace_eq;

      GT* deg_table;

      // To be removed
      int workspace_size;

      size_t small_mgs_size;

      // GT* all;
      // int n;
      // size_t size_all;

      GPUWorkspace
       ( int mon_pos_size, int n_coef, int n_constant, int n_eq, int dim,
         int n_predictor, CT alpha=CT(1,0), int base_table_size=0,
         int n_path=1, int verbose = 0 )
      {
         this->mon_pos_size = mon_pos_size;
         this->dim = dim;
         this->n_eq = n_eq;
         this->n_coef = n_coef;
         this->n_constant = n_constant;
         this->n_predictor = n_predictor;
         this->n_path = n_path;
         if(verbose > 0)
         {
            std::cout << "this->n_path = " << this->n_path << std::endl;
         }
         this->n_path_continuous = n_path;
         n_array = n_predictor + 1;

         this->alpha = alpha;
         cudaMalloc((void **)&alpha_gpu, sizeof(GT));
         GT* alpha_gpu_host = (GT *)malloc(sizeof(GT));
         comp1_qd2gqd(&alpha, alpha_gpu_host);
         cudaMemcpy(alpha_gpu,alpha_gpu_host,sizeof(GT),
                    cudaMemcpyHostToDevice);

         n_row = n_eq;
         n_col = dim+1;

         // Matrix size
         n_matrix = n_row*n_col;
         n_matrix_R = (n_col+1)*n_col/2;
         int row_block = 0;
         if(n_row > BS_QR)
         {
            row_block = (n_row-1)/matrix_block_row+1;
         }
         n_matrix_P = row_block*matrix_block_pivot_col*n_col;

         // GT arrays
         int n_eval_arrays = n_coef+mon_pos_size+n_eq*(dim+1);
         int n_qr_arrays = n_matrix + n_matrix_R + n_matrix_P + dim;
         int n_predict_arrays = dim*(n_predictor+1)+(n_predictor+1);
         int n_x_t_arrays = dim+5;

         n_GT_arrays
          = n_path*(n_eval_arrays+n_qr_arrays+n_predict_arrays+n_x_t_arrays);
         size_GT_arrays = n_GT_arrays*sizeof(GT);
         if(verbose > 0)
         {
            std::cout << "size_GT_arrays = " << size_GT_arrays << std::endl;
         }
         cudaMalloc((void **) &GT_arrays, size_GT_arrays);

         // GT arrays: Eval arrays
         coef   = GT_arrays;
         mon    = coef + n_coef*n_path;
         sum    = mon - n_constant*n_path;
         matrix = mon + mon_pos_size*n_path;
         f_val  = matrix + n_eq*dim*n_path;

         // GT arrays: QR arrays
         if(verbose > 0)
         {
            std::cout << "n_matrix = " << n_matrix
                      << ", n_matrix_R = " << n_matrix_R << std::endl;
         }
         V = matrix;
         R = V+n_matrix*n_path;
         P = R + n_matrix_R*n_path;
         sol = P+n_matrix_P*n_path;
         matrix_horizontal_mult = sol + dim*n_path;

         // GT arrays: Predict arrays
         x_array = matrix_horizontal_mult + n_matrix*n_path;
         t_array = x_array + dim*(n_predictor+1)*n_path;

         // GT arrays: x and t arrays
         one_minor_t   = t_array + (n_predictor+1)*n_path;
         // only for multiple path
         x_mult        = one_minor_t + n_path;
         t_mult        = x_mult + dim*n_path;
         t_last_mult   = t_mult + n_path;
         delta_t_mult  = t_mult + 2*n_path;
         newton_t_mult = t_mult + 3*n_path;

         x = x_array;
         t = t_array;
         x_last = x;
         t_last = t;
         x_t_idx = 0;

         // int arrays
         cudaMalloc((void **) &int_arrays, 7*n_path*sizeof(int));
         x_t_idx_mult   = int_arrays;
         n_point_mult   = int_arrays + n_path;
         path_success   = int_arrays + 2*n_path;
         n_success      = int_arrays + 3*n_path;
         newton_success = int_arrays + 4*n_path;
         path_idx       = int_arrays + 5*n_path;
         end_range      = int_arrays + 6*n_path;

         // double arrays
         cudaMalloc((void **) &double_arrays, 6*n_path*sizeof(double));
         max_delta_x_gpu    = double_arrays;
         r_max_delta_x_gpu  = double_arrays + n_path;
         max_f_val_gpu      = double_arrays + 2*n_path;
         max_f_val_last_gpu = double_arrays + 3*n_path;
         r_max_f_val_gpu    = double_arrays + 4*n_path;
         max_x_gpu          = double_arrays + 5*n_path;

         // host arrays
         path_success_host   = (int *)malloc(n_path*sizeof(int));
         newton_success_host = (int *)malloc(n_path*sizeof(int));
         path_idx_host       = (int *)malloc(n_path*sizeof(int));

         r_max_f_val_host   = (double *)malloc(n_path*sizeof(double));
         max_f_val_host     = (double *)malloc(n_path*sizeof(double));
         max_delta_x_host   = (double *)malloc(n_path*sizeof(double));
         r_max_delta_x_host = (double *)malloc(n_path*sizeof(double));
         max_x_host         = (double *)malloc(n_path*sizeof(double));

         // eq in one block, developping
         workspace_eq=NULL;
         if(verbose > 0)
         {
            std::cout << "GPU initialized" << std::endl;
         }
         // To be removed
         // n = mon_pos_size+n_coef;
         if(verbose > 0)
         {
            std::cout << " mon_pos_size = " << mon_pos_size
                      << " n_coef = " << n_coef << std::endl;
         }
         // size_all = n*sizeof(GT);
         // all = NULL;
         // cudaMalloc((void **)&all, size_all);

         workspace_size = n_matrix;
         workspace_size += n_matrix_R;
         workspace_size += dim;

         small_mgs_size = (n_matrix+n_matrix_R+32)*sizeof(GT) + dim*sizeof(T);

         deg_table = NULL;
         if(verbose > 0)
         {
            std::cout << "base_table_size = " << base_table_size << std::endl;
         }
         if(base_table_size > 0)
         {
            cudaMalloc((void **)&deg_table, base_table_size*sizeof(GT));
         }
      }

      ~GPUWorkspace()
      {
         // std::cout << "Delete GPUWorkspace" << std::endl;
         cudaFree(int_arrays);
      }

      void init_workspace_eq ( int n_pos_total_eq, int n_path );

      void init_matrix ( int dim, int n_eq );

      void init_V_value ( CT* V_cpu );

      void init_x_t ( int dim, int n_predictor );

      void update_x_t_idx();

      void update_x_t_idx_all ( int* x_t_idx_host );

      int* get_x_t_idx_all();

      void update_t_value ( CT cpu_t );

      void update_t_value_array ( CT* cpu_t, int* x_t_idx_host );

      void update_t_value_mult ( CT* cpu_t );

      void update_t_value_inverse ( CT cpu_one_minor_t );

      void update_x_value ( CT* cpu_sol0);

      void update_x_value_array ( CT* cpu_sol0 );

      void update_x_value_mult ( CT* cpu_sol0 );

      void update_x_mult_horizontal ( CT* cpu_sol0 );

      void update_x_t_value ( CT* cpu_sol0, CT cpu_t );

      void update_x_t_value_array
              ( CT* cpu_sol0, CT* cpu_t, int* x_t_idx_host );

      void update_x_t_value_mult ( CT* cpu_sol0, CT* cpu_t );

      void update_x_mult_vertical ( CT* cpu_sol0, int* x_t_idx_host );

      void update_t_value_mult2 ( CT* cpu_t, int* x_t_idx_host );

      void update_x_t ( CT* cpu_sol0, CT cpu_t );

      CT* get_matrix();

      CT* get_matrix ( int sys_idx );

      CT** get_matrix_mult();

      CT* get_workspace ( int sys_idx );

      CT* get_workspace();

      CT** get_workspace_mult();

      CT* get_matrix_r ( int sys_idx=0 );

      void print_matrix_r();

      CT* get_x();

      CT** get_mult_x_horizontal();

      CT* get_x_array();

      CT* get_t_array();

      CT** get_x_all();

      CT** get_x_last_all();

      void print_x();

      void print_x_mult ( int path_idx_one=-1 );

      void print_t_mult ( int path_idx_one=-1 );

      void print_delta_t_mult ( int path_idx_one=-1 );

      void print_x_last_mult();

      void print_x_array();

      void print_t_array();

      CT* get_f_val();

      void print_f_val();

      CT* get_coef_mult();

      CT* get_mon_mult();

      CT* get_x_last();

      CT* get_sol ( int path_idx=0 );

      CT* get_sol_array();

      CT** get_sol_mult();

      T1 sol_norm();

      void init_x_t_predict_test();
};

class GPUInst
{
   public:

      bool PED_hom;
      int n_path;

      // Sol Instruction
      int dim;
      int n_eq;

      /**** workspace Instruction ***/
      int n_workspace;
      int n_constant;

      // Coef Instruction
      int n_coef;
      GT* coef;

      int coef_BS;
      dim3 coef_grid;

      int dim_BS;
      dim3 dim_grid;

      /**** Mon Instruction ****/
      // for leveled kernels
      int level;
      int* n_mon_level;
      // for single kernel
      int n_mon;
      int* mon_pos_start;
      unsigned short* mon_pos;

      int n_mon_global;

      dim3* mon_level_grid;
      int* n_mon_level_rest;
      dim3* mon_level_grid_rest;

      int mon_pos_size;

      int mon_level0_BS;
      int mon_level_BS;

      int mon_global_BS;
      dim3 mon_global_grid;

      int n_mon_block;
      dim3 mon_block_grid;
      int BS_mon_block;
      int NB_mon_block;
      int* mon_pos_start_block;
      unsigned short* mon_pos_block;
      int n_mon_single;
      unsigned short* mon_single_pos_block;

      /**** Sum instruction ****/
      int n_sum; // size of sum_start
      int n_sum_levels;
      int* n_sum_level;
      dim3* sum_level_grid;
      int* n_sum_level_rest;
      dim3* sum_level_grid_rest;

      int* sum_pos_start;
      int* sum_pos;

      int* sum_pos_start_align;
      int* sum_pos_align;

      int sum_BS;
      dim3 sum_grid;

      int n_step_GPU;
      int n_point_GPU;
      int n_eval_GPU;
      int n_mgs_GPU;

      int predict_BS;
      dim3 predict_grid;

      int* eq_pos_start;
      int n_mon_total_eq;
      int* mon_pos_start_eq;
      GT* coef_eq;
      int n_pos_total_eq;
      unsigned short* mon_pos_eq;

      CT alpha;

      int n_sum_zero;
      int* sum_zeros;

      int base_table_size;
      int* base_table_start;
      int* max_deg_base;
      unsigned short* mon_exp;
      int n_mon_base;
      int n_mon_base_start;

      GPUInst ( const CPUInstHom& cpu_inst, int n_path )
      {
         PED_hom = cpu_inst.PED_hom;
         dim = cpu_inst.dim;
         n_eq = cpu_inst.n_eq;
         this->n_path = n_path;
         init_predict();
         init_coef(cpu_inst.CPU_inst_hom_coef);
         if(MON_EVAL_METHOD == 1 && n_path == 1)
         {
            init_mon(cpu_inst.CPU_inst_hom_block);
            init_sum(cpu_inst.CPU_inst_hom_sum_block,
                     cpu_inst.CPU_inst_hom_sum);
         }
         else
         {
            init_mon(cpu_inst.CPU_inst_hom_mon);
            init_sum(cpu_inst.CPU_inst_hom_sum);
         }
         init_workspace(cpu_inst);

         dim_BS = 32;
         dim_grid = get_grid(dim,dim_BS,n_path);
         n_step_GPU = 0;
         n_point_GPU = 0;
         n_eval_GPU = 0;
         n_mgs_GPU = 0;

         init_eq(cpu_inst.CPU_inst_hom_eq);

         // Initialize the base part
         if(cpu_inst.CPU_inst_hom_mon.max_deg_base != NULL)
         {
            init_base(cpu_inst.CPU_inst_hom_mon);
	 }
	 else
         {
            // No base
            base_table_size = 0;
            base_table_start = NULL;
            max_deg_base = NULL;
            mon_exp = NULL;
            n_mon_base = 0;
            n_mon_base_start = 0;
         }
      }

      ~GPUInst()
      {
          cudaFree(coef);
          cudaFree(mon_pos_start);
          cudaFree(mon_pos);
          cudaFree(sum_pos_start);
          cudaFree(sum_pos);
      }

      void init_predict();

      void init_coef ( const CPUInstHomCoef& cpu_inst_coef );

      void init_mon ( const CPUInstHomMon& cpu_inst_mon );

      void init_mon ( const CPUInstHomMonBlock& cpu_inst_mon_block );

      void init_base ( const CPUInstHomMon& cpu_inst_mon);

      void init_sum ( const CPUInstHomSumBlock& cpu_inst_sum,
                      const CPUInstHomSum& cpu_inst_sum_orig );

      void init_sum ( const CPUInstHomSum& cpu_inst_sum );

      void init_workspace ( const CPUInstHom& cpu_inst );

      void init_eq ( const CPUInstHomEq& cpu_inst_mon_eq );
};

#endif /* GPU_DATA_H_ */
