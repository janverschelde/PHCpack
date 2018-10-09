// The file cpuinsthomsumblock.h defines the class CPUInstHomSumBlock.

#ifndef __CPUINSTHOMSUMBLOCK_H__
#define __CPUINSTHOMSUMBLOCK_H__

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

#include "cpuinsthomsumblock.tpp"

#endif /* __CPUINSTHOMMON_H__ */
