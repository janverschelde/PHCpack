// The file cpuinsthommon.h defines the class CPUInstHomMon.

#ifndef __CPUINSTHOMMON_H__
#define __CPUINSTHOMMON_H__

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

#include "cpuinsthommon.tpp"

#endif /* __CPUINSTHOMMON_H__ */
