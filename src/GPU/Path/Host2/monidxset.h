// The file monidxset.h defines the class MonIdxSet.

#ifndef __MONIDXSET_H__
#define __MONIDXSET_H__

#include<iostream>
#include<algorithm>

#include "twoint.h"
#include "generalset.h"

template <class ComplexType>
class MonIdxSet
{
   // friend class MonSet;
   template <class T> friend class MonSet;

   IntSet pos;
   IntSet exp;
   int eq_idx;
   int mon_idx;
   bool sys_idx;
   ComplexType coef;

   void init ( const int n, const int* pos, const int* exp )
   {
      this->pos = IntSet(n, pos);
      this->exp = IntSet(n, exp);
      eq_idx = 0;
      mon_idx = 0;
      sys_idx = 0;
      coef = 0.0;
   }

   void init ( const int n, const int* pos, const int* exp,
               int eq_idx, int mon_idx, bool sys_idx )
   {
      this->pos = IntSet(n, pos);
      this->exp = IntSet(n, exp);
      this->eq_idx = eq_idx;
      this->mon_idx = mon_idx;
      this->sys_idx = sys_idx;
      this->coef = 0.0;
   }

   void init ( const int n, const int* pos, const int* exp,
               int eq_idx, int mon_idx, bool sys_idx,
               const ComplexType& coef )
   {
      this->pos = IntSet(n, pos);
      this->exp = IntSet(n, exp);
      this->eq_idx = eq_idx;
      this->mon_idx = mon_idx;
      this->sys_idx = sys_idx;
      this->coef = coef;
   }

   void init ( const IntSet& pos, const IntSet& exp,
               int eq_idx, int mon_idx, bool sys_idx,
               const ComplexType& coef )
   {
      this->pos = pos;
      this->exp = exp;
      this->eq_idx = eq_idx;
      this->mon_idx = mon_idx;
      this->sys_idx = sys_idx;
      this->coef = coef;
   }

   public:

      MonIdxSet()
      {
         eq_idx = 0;
         mon_idx = 0;
         sys_idx = 0;
         coef = 0.0;
      };

      MonIdxSet ( const int n, const int* pos, const int* exp )
      {
         init(n, pos, exp);
      }

      MonIdxSet ( const int n, const int* pos, const int* exp,
                 int eq_idx, int mon_idx, bool sys_idx )
      {
         init(n, pos, exp, eq_idx, mon_idx, sys_idx);
      }

      MonIdxSet ( const int n, const int* pos, const int* exp,
                  int eq_idx, int mon_idx, bool sys_idx,
                  const ComplexType& coef )
      {
         init(n, pos, exp, eq_idx, mon_idx, sys_idx, coef);
      }

      MonIdxSet ( const MonIdxSet& original )
      {
         init(original.pos, original.exp, original.eq_idx,
              original.mon_idx, original.sys_idx, original.coef);
      }

      MonIdxSet& operator= ( const MonIdxSet& original )
      {
         init(original.pos, original.exp, original.eq_idx,
              original.mon_idx, original.sys_idx, original.coef);
         return *this;
      }

      ~MonIdxSet()
      {
         //std::cout << "delete MonIdxSet " << this << endl;
      }

      // Basic functions

      IntSet get_pos()
      {
         return pos;
      }

      IntSet get_exp()
      {
         return exp;
      }

      int get_eq_idx()
      {
         return eq_idx;
      }

      int get_mon_idx()
      {
         return mon_idx;
      }

      bool get_sys_idx()
      {
         return sys_idx;
      }

      ComplexType get_coef()
      {
         return coef;
      }

      // comparison operators

      bool operator < ( const MonIdxSet& that ) const
      {
         // Compare exponent set
         if(this->exp < that.exp)
            return true;
         else if(this->exp > that.exp)
            return false;

         // Compare exponent set
         if(this->pos < that.pos)
            return true;
         else if(this->pos > that.pos)
            return false;

         // Compare equation index
         if(this->eq_idx < that.eq_idx)
            return true;
         else if(this->eq_idx > that.eq_idx)
            return false;

         // Compare system index
         if(this->sys_idx < that.sys_idx)
            return true;
         else
            return false;
      }

      bool operator > ( const MonIdxSet& that ) const
      {
         // Compare exponent set
         if(this->exp > that.exp)
            return true;
         else if(this->exp < that.exp)
            return false;

         // Compare exponent set
         if(this->pos > that.pos)
            return true;
         else if(this->pos < that.pos)
            return false;
         //
         // Compare equation index
         if(this->eq_idx > that.eq_idx)
            return true;
         else if(this->eq_idx < that.eq_idx)
            return false;

         // Compare system index
         if(this->sys_idx > that.sys_idx)
            return true;
         else
            return false;
      }

      bool operator == ( const MonIdxSet& that ) const
      {
         if(this->pos == that.pos && this->exp == that.exp)
            return true;
         else
            return false;
      }

      void print()
      {
         pos.print();
         exp.print();
         std::cout << "s" << sys_idx
                   << " e" << eq_idx
                   << " m" << mon_idx << std::endl;
         std::cout << coef;
         std::cout << std::endl;
      }

      void sorted()
      {
         int n_elements =pos.get_n();
         TwoInt* tmp_pos_exp = new TwoInt[n_elements];

         for(int idx=0; idx<n_elements; idx++)
         {
            tmp_pos_exp[idx].init(pos.get_element(idx), exp.get_element(idx));
         }
         PosExpSet tmp_set(n_elements, tmp_pos_exp);
         tmp_set.sorted();

         int* tmp_pos = new int[n_elements];
         int* tmp_exp = new int[n_elements];
         for(int idx=0; idx<n_elements; idx++)
         {
            TwoInt tmpTwoInt = tmp_set.get_element(idx);
            tmp_pos[idx] = tmpTwoInt.int1;
            tmp_exp[idx] = tmpTwoInt.int2;
         }
         pos.update(n_elements, tmp_pos);
         exp.update(n_elements, tmp_exp);
      }
};

#endif /* __MONIDXSET_H__ */
