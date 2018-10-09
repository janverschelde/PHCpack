// The file monset.h defines the class MonSet.

#ifndef __MONSET_H__
#define __MONSET_H__

#include<iostream>
#include<algorithm>

#include "monidxset.h"
#include "eqidxcoef.h"

template <class ComplexType>
class MonSet
{
   IntSet pos;
   IntSet exp;
   EqSet eq_idx;

   void init ( const int n, const int* pos, const int* exp )
   {
      this->pos = IntSet(n, pos);
      this->exp = IntSet(n, exp);
   }

   void init ( const MonSet& original )
   {
      pos = original.pos;
      exp = original.exp;
      eq_idx = original.eq_idx;
   }

   public:

      MonSet(){ };

      void copy_pos ( const MonIdxSet<ComplexType>& original )
      {
         this->pos = original.pos;
         this->exp = original.exp;
      }

      void update_eq_idx ( int n, EqIdxCoef<ComplexType>* eq_idx_elements )
      {
         eq_idx.update(n,eq_idx_elements);
      }

      MonSet ( const int n, const int* pos, const int* exp )
      {
         init(n, pos, exp);
      }

      MonSet ( const MonSet& original )
      {
         init(original);
      }

      MonSet& operator= ( const MonSet& original )
      {
         init(original);
         return *this;
      }

      ~MonSet(){ }

      bool operator < ( const MonSet& that ) const
      {
         if(this->exp < that.exp)
            return true;
         else if(this->exp > that.exp)
            return false;
         if(this->pos < that.pos)
            return true;
         else
            return false;
      }

      bool operator > ( const MonSet& that ) const
      {
         if(this->exp > that.exp)
            return true;
         else if(this->exp < that.exp)
            return false;
         if(this->pos > that.pos)
            return true;
         else
            return false;
      }

      bool operator == ( const MonSet& that ) const
      {
         if(this->pos == that.pos && this->exp == that.exp)
            return true;
         else
            return false;
      }

      void sorted()
      {
         int n_elements =pos.get_n();
         TwoInt* tmp_pos_exp = new TwoInt[n_elements];

         for(int idx=0; idx<n_elements; idx++)
         {
            tmp_pos_exp[idx].init(pos.get_element(idx),
                                  exp.get_element(idx));
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

      int get_n()
      {
         return pos.get_n();
      }

      int get_n_mon()
      {
         return eq_idx.get_n();
      }

      int get_pos ( int pos_idx )
      {
         return pos.get_element(pos_idx);
      }

      int get_exp ( int pos_idx )
      {
         return exp.get_element(pos_idx);
      }

      int get_eq_idx ( int idx )
      {
         return eq_idx.elements[idx].get_eq_idx();
      }

      void write_coef ( ComplexType*& tmp_coef )
      {
         for(int i=0; i<eq_idx.n; i++)
            eq_idx.elements[i].write_coef(tmp_coef);
      }

      void write_coef ( ComplexType*& tmp_coef, int mon_idx )
      {
         eq_idx.elements[mon_idx].write_coef(tmp_coef);
      }

      void write_pos ( int*& tmp_pos )
      {
         for(int i=0; i<pos.n; i++) (*tmp_pos++) = pos.elements[i];
      }

      void write_pos ( unsigned short*& tmp_pos )
      {
         for(int i=0; i<pos.n; i++) (*tmp_pos++) = pos.elements[i];
      }

      void write_exp ( int*& tmp_exp )
      {
         for(int i=0; i<exp.n; i++) (*tmp_exp++) = exp.elements[i];
      }

      void write_exp ( unsigned short*& tmp_exp )
      {
         for(int i=0; i<exp.n; i++) (*tmp_exp++) = exp.elements[i];
      }

      bool has_base()
      {
         for(int i=0; i<exp.n; i++)
            if(exp.elements[i]>1) return true;

         return false;
      }

      int get_base_size()
      {
         for(int i=0; i<exp.n; i++)
            if(exp.elements[i] > 1) return exp.n - i;

         return 0;
      }

      int get_base_start()
      {
         for(int i=0; i<exp.n; i++)
            if(exp.elements[i] > 1) return i;

         return exp.n;
      }

      friend std::ostream& operator << ( std::ostream& o, const MonSet& c )
      {
         return o << c.pos << c.exp << c.eq_idx << std::endl;
      }
};

#endif /* __MONSET_H__ */
