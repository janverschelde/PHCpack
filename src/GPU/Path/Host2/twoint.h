// The file twoint.h defines the class TwoInt to store a pair of integers.

#ifndef __TWOINT_H__
#define __TWOINT_H__

class TwoInt
{
   public:

      int int1;
      int int2;

      TwoInt()
      {
         int1 = 0;
         int2 = 0;
      }

      TwoInt ( const int int1, const int int2 )
      {
         init(int1, int2);
      }

      void init ( int int1, int int2 )
      {
         this->int1 = int1;
         this->int2 = int2;
      }

      bool operator < ( const TwoInt& that ) const
      {
         if(this->int2 < that.int2)
            return true;
         else if(this->int2 > that.int2)
            return false;

         if(this->int1 < that.int1)
            return true;
         else
            return false;
      }

      bool operator > ( const TwoInt& that ) const
      {
         if(this->int2 > that.int2)
            return true;
         else if(this->int2 < that.int2)
            return false;

         if(this->int1 > that.int1)
            return true;
         else
            return false;
      }

      bool operator == ( const TwoInt& that ) const
      {
         if(this->int1 == that.int1 && this->int2 == that.int2)
            return true;
         else
            return false;
      }
};

#endif /* __TWOINT_H__ */
