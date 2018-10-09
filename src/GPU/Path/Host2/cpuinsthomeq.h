// The file cpuinsthomeq defines the class CPUInstHomEq.

#ifndef __CPUINSTHOMEQ_H__
#define __CPUINSTHOMEQ_H__

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

#include "cpuinsthomeq.tpp"

#endif /* __CPUINSTHOMEQ_H__ */
