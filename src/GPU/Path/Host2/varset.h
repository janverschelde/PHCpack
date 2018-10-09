// The file varset.h defines the class VarSet to store elements.

#ifndef __VARSET_H__
#define __VARSET_H__

#include<iostream>
#include<algorithm>

class VarSet
{
   public:

      int n;
      VarSet* elements;

      void init()
      {
         this->n = 0;
         elements = NULL;
      }

      void init ( const int n )
      {
         this->n = n;
         elements = NULL;
      }

      void init ( const int n, const VarSet* elements )
      {
         this->n = n;
         if(elements != NULL)
         {
            this->elements = new VarSet[n];
            for(int i=0; i<n; i++) this->elements[i] = elements[i];
         }
         else
         {
            this->elements = NULL;
         }
      };

      VarSet()
      {
         init();
      }

      VarSet( int n )
      {
         init(n);
      }

      VarSet ( const int n, const VarSet* elements )
      {
         init(n, elements);
      }

      VarSet ( const VarSet& original )
      {
         init(original.n, original.elements);
      }

      VarSet& operator= ( const VarSet& original )
      {
         init(original.n, original.elements);
         return *this;
      }

      bool operator < ( const VarSet& that ) const
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
            if(this->elements == NULL && that.elements == NULL)
            {
               return false;
            }
            else if(this->elements == NULL)
            {
               std::cout << "Not on the same level" << std::endl;
               return false;
            }
            else if(that.elements == NULL)
            {
               std::cout << "Not on the same level" << std::endl;
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

      bool operator > ( const VarSet& that ) const
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
            if(this->elements == NULL && that.elements == NULL)
            {
               return false;
            }
            else if(this->elements == NULL)
            {
               std::cout << "Not on the same level" << std::endl;
               return false;
            }
            else if(that.elements == NULL)
            {
               std::cout << "Not on the same level" << std::endl;
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

      ~VarSet()
      {
         delete[] elements;
         elements = NULL;
      }

      void print ( int level = 0 )
      {
         if(elements != NULL)
         {
            std::cout << std::endl;
            for(int i=0; i<level; i++) std::cout << "   ";

            std::cout << "l"<< level << " n" << n << " : ";
            for(int i=0; i<n; i++) elements[i].print(level+1);
         }
         else
         {
            std::cout << n << ", ";
         }
         if(level == 0) std::cout << std::endl;
      }

      void sorted()
      {
         if(elements != NULL)
         {
            for(int i=0; i<n; i++) elements[i].sorted();

            std::sort(elements, elements+n);
         }
      }
};

#endif /* __VARSET_H__ */
