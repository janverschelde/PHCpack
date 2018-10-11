// The file cpuinsthom.h defines the class CPUInstHom.

#ifndef __CPUINSTHOM_H__
#define __CPUINSTHOM_H__

template <class ComplexType, class RealType>
class CPUInstHom
{
   public:

      bool PED_hom; // true: homotopy, false: single
      CPUInstHomCoef<ComplexType,RealType> CPU_inst_hom_coef;
      CPUInstHomMon<ComplexType> CPU_inst_hom_mon;
      CPUInstHomSum<ComplexType> CPU_inst_hom_sum;
      CPUInstHomMonBlock<ComplexType> CPU_inst_hom_block;
      CPUInstHomSumBlock<ComplexType> CPU_inst_hom_sum_block;

      CPUInstHomEq<ComplexType> CPU_inst_hom_eq;

      int n_constant;
      int dim;
      int n_eq;
      int n_coef;
      int n_predictor;

      // Record timing for both CPU and GPU
      Path<ComplexType,RealType> path_data;
      Path<ComplexType,RealType> path_data_gpu;
      double timeSec_Path_CPU;
      double timeSec_Path_GPU;

      bool success_CPU;
      bool success_GPU;
      ComplexType t_CPU;
      ComplexType t_GPU;
      RealType max_residual;
      RealType max_delta_x;

      int n_step_CPU;
      int n_step_GPU;
      int n_point_CPU;
      int n_point_GPU;
      int n_eval_CPU;
      int n_eval_GPU;
      int n_mgs_CPU;
      int n_mgs_GPU;

      void init ( MonSet<ComplexType>* hom_monset, int n_monset,
                  int n_constant, int total_n_mon, int dim,
                  int n_eq, int n_predictor, ComplexType alpha,
                  int* max_deg_base, int verbose = 0 );

      void init ( PolySys<ComplexType,RealType>& Target_Sys,
                  PolySys<ComplexType,RealType>& Start_Sys, int dim, int n_eq,
                  int n_predictor, ComplexType alpha, int verbose = 0 );

      void init ( PolySys<ComplexType,RealType>& Target_Sys,
                  int dim, int n_eq, int n_predictor,
                  ComplexType alpha, int verbose = 0 );

      CPUInstHom()
      {
         PED_hom = false;
         n_constant = 0;
         dim = 0;
         n_eq = 0;
         n_coef = 0;
         n_predictor = 0;
         // For record only
         timeSec_Path_CPU = 0;
         timeSec_Path_GPU = 0;
         success_CPU = 0;
         success_GPU = 0;
         t_CPU = ComplexType(0.0,0.0);
         t_GPU = ComplexType(0.0,0.0);
         max_residual = 0.0;
         max_delta_x = 0.0;
         n_step_CPU = 0;
         n_step_GPU = 0;
         n_point_CPU = 0;
         n_point_GPU = 0;
         n_eval_CPU = 0;
         n_eval_GPU = 0;
         n_mgs_CPU = 0;
         n_mgs_GPU = 0;
      }

      CPUInstHom
       ( MonSet<ComplexType>* hom_monset, int n_monset, int n_constant,
         int total_n_mon, int dim, int n_eq, int n_predictor,
         ComplexType alpha, int* max_deg_base )
      {
         // For record only
         timeSec_Path_CPU = 0;
         timeSec_Path_GPU = 0;
         success_CPU = 0;
         success_GPU = 0;
         t_CPU = ComplexType(0.0,0.0);
         t_GPU = ComplexType(0.0,0.0);
         max_residual = 0.0;
         max_delta_x = 0.0;
         n_step_CPU = 0;
         n_step_GPU = 0;
         n_point_CPU = 0;
         n_point_GPU = 0;
         n_eval_CPU = 0;
         n_eval_GPU = 0;
         n_mgs_CPU = 0;
         n_mgs_GPU = 0;
         init(hom_monset,n_monset,n_constant,total_n_mon,
              dim,n_eq,n_predictor,alpha,max_deg_base);
      }

      ~CPUInstHom(){ }

      void print();
 
      void init_workspace
       ( Workspace<ComplexType>& workspace_cpu, int verbose=0 );
      /*
         Initialization of the data structures to evaluate and
         differentiate the polynomial system.
         This function must be executed before the eval.
       */

      void eval ( Workspace<ComplexType>& workspace_cpu,
                  const ComplexType* sol,
                  const ComplexType t, int reverse=0 );
      /*
         Evaluates the polynomials corresponding to the workspace
         at the array of complex numbers in sol, at t if a homotopy.
         The number of elements in sol should equal the dimension
         of the polynomial system.
         The result is in workspace_cpu.matrix, which is a single
         index array of (dim+1)*dim complex numbers, where dim is
         the dimension of the polynomial system.
         The function value of the system at sol is in the last dim
         elements of the matrix.  The first dim*dim numbers contain
         the values of the Jacobian matrix at sol.
       */
      
      void update_alpha ( ComplexType alpha=ComplexType(0.0,0.0) );

      void compare_path_cpu_gpu();

};

#include "cpuinsthom.tpp"

#endif /* __CPUINSTHOM_H__ */
