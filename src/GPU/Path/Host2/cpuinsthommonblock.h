// The file cpuinsthommonblock.h defines the class CPUInstHomMonBlock.

#ifndef __CPUINSTHOMMONBLOCK_H__
#define __CPUINSTHOMMONBLOCK_H__

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

#include "cpuinsthommonblock.tpp"

#endif /* __CPUINSTHOMMONBLOCK_H__ */
