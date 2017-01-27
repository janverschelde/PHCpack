// eval_host.h contains templated data structure for the evaluation and
// differentiation of polynomials on the host.
// It includes the file eval_host.tpp.

#ifndef CPU_INSTRUCTION_EVAL_H_
#define CPU_INSTRUCTION_EVAL_H_

#include <sys/time.h>
#include <unistd.h>
#include "varset.h"
#include "parameter.h"
#include "workspace_host.h"
#include "path_data.h"
#include "poly.h"

#define warp_size 32

template <class ComplexType, class RealType>
class CPUInstHomCoef
{
   public:

      int n_coef;
      ComplexType* coef_orig;
      ComplexType alpha;

      CPUInstHomCoef()
      {
         n_coef = 0;
         coef_orig = NULL;
         alpha = ComplexType(0.0,0);
      }

      CPUInstHomCoef ( MonSet<ComplexType>* hom_monset,
                       int total_n_mon, int n_monset,
                       int n_constant, ComplexType alpha )
      {
         init(hom_monset, total_n_mon, n_monset, n_constant, alpha);
      }

      ~CPUInstHomCoef()
      {
         // std::cout << "Delete CPUInstHomCoef" << std::endl;
         if(coef_orig != NULL) delete[] coef_orig;
      }

      void init ( MonSet<ComplexType>* hom_monset,
                  int total_n_mon, int n_monset,
                  int n_constant, ComplexType alpha, int verbose = 0 );

      void print();

      void eval ( const ComplexType t, ComplexType* coef, int reverse=0 );

      void update_alpha ( ComplexType alpha=ComplexType(0.0,0.0) );
};

template <class ComplexType>
class CPUInstHomMon
{
   public:

      int level; // size of n_mon_level
      int* n_mon_level;
      int n_mon; // size of pos_start, sum of n_mon_level
      int* mon_pos_start;
      int mon_pos_size;
      unsigned short* mon_pos;
      unsigned short* mon_exp;
      int* max_deg_base;
      int n_mon_base_start;
      int flops_multiple;

      CPUInstHomMon()
      {
         level = 0;
         n_mon_level = NULL;
         n_mon = 0;
         mon_pos_start = NULL;
         mon_pos_size = 0;
         mon_pos = NULL;
         mon_exp = NULL;
         max_deg_base = NULL;
         flops_multiple = 0;
         n_mon_base_start = 0;
      }

      CPUInstHomMon ( MonSet<ComplexType>* hom_monset,
                      int n_monset, int total_n_mon,
                      int n_constant, int* max_deg_base, int verbose = 0 )
      {
         CPUInstHomMon();
         init(hom_monset,n_monset,total_n_mon,n_constant,max_deg_base,verbose);
      }

      ~CPUInstHomMon()
      {
         // std::cout << "Delete CPUInstHomMon" << std::endl;
         if(n_mon_level != NULL) delete[] n_mon_level;
         if(mon_pos_start != NULL) delete[] mon_pos_start;
         if(mon_pos != NULL) delete[] mon_pos;
         if(mon_exp != NULL) delete[] mon_exp;
         if(max_deg_base != NULL) delete[] max_deg_base;
      }

      void init ( MonSet<ComplexType>* hom_monset, int n_monset,
                  int total_n_mon,
                  int n_constant, int* max_deg_base, int verbose = 0 );

      void eval ( int dim, const ComplexType* x_val, ComplexType* mon,
                  ComplexType* coef, ComplexType** deg_table=NULL );

      void eval_deg_table ( int dim, const ComplexType* x_val,
                            ComplexType** deg_table );

      void eval_base ( ComplexType** deg_table, ComplexType* coef );

      void print();
};

template <class ComplexType>
class CPUInstHomMonBlock
{
   public:

      int n_mon;
      int BS;
      int NB;
      int mon_pos_block_size;
      int* mon_pos_start_block;
      unsigned short* mon_pos_block;
      unsigned short* max_var_block;
      int n_mon_single;
      unsigned short* mon_single_pos_block;

      CPUInstHomMonBlock()
      {
         n_mon = 0;
         n_mon_single = 0;
         BS = 0;
         NB = 0;
         mon_pos_block_size = 0;
         mon_pos_start_block = NULL;
         mon_pos_block = NULL;
         max_var_block = NULL;
         mon_single_pos_block = NULL;
      }

      CPUInstHomMonBlock
       ( CPUInstHomMon<ComplexType>& orig, int BS, int verbose )
      {
         init(orig,BS,verbose);
      }

      ~CPUInstHomMonBlock()
      {
         if(max_var_block != NULL) delete[] max_var_block;
         if(mon_pos_start_block != NULL) delete[] mon_pos_start_block;
         if(mon_pos_block != NULL) delete[] mon_pos_block;
      }

      void init ( CPUInstHomMon<ComplexType>& orig, int BS, int verbose );

      void print();
};

template <class ComplexType>
class CPUInstHomSumBlock
{
   public:

      int n_sum; // size of sum_start
      int n_sum_levels;
      int* n_sum_level;
      int* n_sum_level_rest;
      int* sum_pos_start;
      int sum_pos_size;
      int* sum_pos;

      CPUInstHomSumBlock()
      {
         n_sum = 0;
         n_sum_levels = 0;
         n_sum_level = NULL;
         n_sum_level_rest = NULL;
         sum_pos_start = NULL;
         sum_pos_size = 0;
         sum_pos = NULL;
      }

      CPUInstHomSumBlock
       ( MonSet<ComplexType>* hom_monset, int n_monset, 
         const int* mon_pos_start, int dim,
         int n_eq, int n_constant, int n_mon0, int* mon_pos_start_block,
         int verbose )
      {
         init(hom_monset,n_monset,mon_pos_start,dim,n_eq,n_constant,n_mon0,
              mon_pos_start_block,verbose);
      }

      ~CPUInstHomSumBlock()
      {
         // std::cout << "Delete CPUInstHomSum" << std::endl;
         if(sum_pos_start != NULL) delete[] sum_pos_start;
         if(sum_pos != NULL) delete[] sum_pos;
      }

      void init ( MonSet<ComplexType>* hom_monset,
                  int n_monset, const int* mon_pos_start,
                  int dim, int n_eq, int n_constant, int n_mon0,
                  int* mon_pos_start_block, int verbose = 0 );

      void eval ( ComplexType* sum, ComplexType* matrix );

      void print();
};

template <class ComplexType>
class CPUInstHomSum
{
   public:

      int n_sum; // size of sum_start
      int n_sum_levels;
      int* n_sum_level;
      int* n_sum_level_rest;
      int* sum_pos_start;
      int sum_pos_size;
      int* sum_pos;
      int n_sum_zero;
      int* sum_zeros;

      CPUInstHomSum()
      {
         n_sum = 0;
         n_sum_levels = 0;
         n_sum_level = NULL;
         n_sum_level_rest = NULL;
         sum_pos_start = NULL;
         sum_pos_size = 0;
         sum_pos = NULL;
         n_sum_zero = 0;
         sum_zeros = NULL;
      }

      CPUInstHomSum
       ( MonSet<ComplexType>* hom_monset,
         int n_monset, const int* mon_pos_start,
         int dim, int n_eq, int n_constant )
      {
         init(hom_monset,n_monset,mon_pos_start,dim,n_eq,n_constant);
      }

      ~CPUInstHomSum()
      {
         // std::cout << "Delete CPUInstHomSum" << std::endl;
         if(sum_pos_start != NULL) delete[] sum_pos_start;
         if(sum_pos != NULL) delete[] sum_pos;
      }

      void init ( MonSet<ComplexType>* hom_monset,
                  int n_monset, const int* mon_pos_start,
                  int dim, int n_eq, int n_constant, int verbose = 0 );

      void eval ( ComplexType* sum, ComplexType* matrix );

      void print();
};

template <class ComplexType>
class CPUInstHomEq
{
   public:

      int n_eq;
      int* n_mon_eq;
      int* eq_pos_start;
      int n_mon_total;
      int* mon_pos_start_eq;
      int n_pos_total;
      unsigned short * mon_pos_eq;
      ComplexType* coef;

      CPUInstHomEq()
      {
         n_eq = 0;
         n_mon_eq = NULL;
         eq_pos_start = NULL;
         n_mon_total = 0;
         mon_pos_start_eq = NULL;
         n_pos_total = 0;
         mon_pos_eq = NULL;
         coef = NULL;
      }

      CPUInstHomEq
       ( MonSet<ComplexType>* hom_monset,
         int n_monset, int n_eq, int n_constant )
      {
         init(hom_monset,n_monset,n_eq,n_constant);
      }

      void init ( MonSet<ComplexType>* hom_monset,
                  int n_monset, int n_eq, int n_constant );

      void print();
};

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
 
      void init_workspace ( Workspace<ComplexType>& workspace_cpu );
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

#include "eval_host.tpp"

#endif /* CPU_INSTRUCTION_EVAL_H_ */
