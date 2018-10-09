// The file generalset.h defines the templated class GeneralSet,
// to store a general set of elements.

#ifndef __GENERALSET_H__
#define __GENERALSET_H__

#include<iostream>
#include<algorithm>

#define IntSet GeneralSet<int>
#define ShortSet GeneralSet<short>
#define EqSet GeneralSet<EqIdxCoef<ComplexType> >
#define PosExpSet GeneralSet<TwoInt>

template <class T>
class GeneralSet
{
   template <class ComplexType> friend class MonSet;
   int n;
   T* elements;

   void init ( const int n, const T* elements )
   {
      this->n = n;
      this->elements = new T[n];
      for(int i=0; i<n; i++) this->elements[i] = elements[i];
   }

   public:

      GeneralSet()
      {
         n = 0;
         elements = NULL;
      };

      GeneralSet ( const int n, const T* elements )
      {
         init(n, elements);
      }

      GeneralSet& operator= ( const GeneralSet& original )
      {
         init(original.n, original.elements);
         return *this;
      }

      ~GeneralSet()
      {
         delete[] elements;
         elements = NULL;
      }

      int get_n()
      {
         return n;
      }

      T get_element ( int i )
      {
         return elements[i];
      }

      void update ( int n, T* elements )
      {
         this->n = n;
         this->elements = elements;
      }

      bool operator < ( const GeneralSet& that ) const
      {
         if(this->n < that.n)
         {
            return true;
         }
         else if(this->n > that.n)
         {
            return false;
         }
         else if(this->n > 0)
         {
            // Check if one is empty
            if(this->elements == NULL || that.elements == NULL)
            {
               std::cout << "At least, one set is empty" << std::endl;
               return false;
            }
            // Compare each element
            for(int i=0; i<n; i++)
            {
               if(this->elements[i] < that.elements[i])
               {
                  return true;
               }
               else if(this->elements[i] > that.elements[i])
               {
                  return false;
               }
            }
            return false;
         }
         else
         {
            return false;
         }
      }

      bool operator > ( const GeneralSet& that ) const
      {
         if(this->n > that.n)
         {
            return true;
         }
         else if(this->n < that.n)
         {
            return false;
         }
         else if(this->n > 0)
         {
            // Check if one is empty
            if(this->elements == NULL || that.elements == NULL)
            {
               std::cout << "At least, one set is empty" << std::endl;
               return false;
            }
            for(int i=0; i<n; i++)
            {
               if(this->elements[i] > that.elements[i])
               {
                  return true;
               }
               else if(this->elements[i] < that.elements[i])
               {
                  return false;
               }
            }
            return false;
         }
         else
         {
            return false;
         }
      }

      bool operator == ( const GeneralSet& that ) const
      {
         if(this->n == that.n)
         {
            for(int i=0; i<n; i++)
               if(this->elements[i] != that.elements[i]) return false;
         }
         else
         {
            return false;
         }
         return true;
      }

      void print()
      {
         std::cout << "n = " << n << ": ";
         for(int i=0; i<n; i++)
            std::cout << elements[i] << "  ";

         std::cout << std::endl;
      }

      void sorted()
      {
         std::sort(elements, elements+n);
      }

      friend std::ostream& operator << ( std::ostream& o,
                                         const GeneralSet<T>& c )
      {
         o << "n = " << c.n << ": ";
         for(int i=0; i<c.n; i++) o << c.elements[i] << "  ";
         o << std::endl;
         return o;
      }
};

#endif /* __GENERALSET_H__ */
