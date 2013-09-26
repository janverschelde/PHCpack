with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Multprec_Natural64_Coefficients is

-- DESCRIPTION :
--   This package provides manipulations on coefficient vectors of natural
--   numbers to support multi-precision arithmetic, based on 64-bit arithmetic.
--   The coefficient vector c represents a multi-precision natural number
--   as the sum of c(i)*Base^i, with c(i) < Base, for i in 0..c'last.
--   The operations in this package are provided for efficiency.
--   No automatic tests are built in to check the requirements.
--   Routines that call these operations must make sure the requirements
--   for the operations are met.

-- DATASTRUCTURE :

  type Array_of_Naturals is array ( natural32 range <> ) of natural64;

-- SELECTORS :

  function Radix return natural64;

  -- DESCRIPTION :
  --   Returns the radix of the number representation.

  function Base return natural64;

  -- DESCRIPTION :
  --   Returns the base = Radix^Exponent of the number representation.

  function Exponent return natural32;

  -- DESCRIPTION :
  --   Returns the exponent which determines the size of the base.

  function Size_of_Coefficient ( n : natural64 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of digits for a coefficient n < Base
  --   or returns Exponent+1 if n is not a coefficient.

  function Digits_to_Normal ( n : Array_of_Naturals ) return natural32;

  -- DESCRIPTION :
  --   An array n is normalized if n(n'last) >= Base/10.
  --   Returns the number of digits need to normalize n.
  --   In case n is entirely zero, zero is returned.

-- COMPARISONS :

  function Equal ( n1,n2 : Array_of_Naturals ) return boolean;

  -- DESCRIPTION :
  --   Returns true if n1 and n2 have the same dimensions and content.

  function "<" ( n1,n2 : Array_of_Naturals ) return boolean;
  function ">" ( n1,n2 : Array_of_Naturals ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the operations hold for their meaning as
  --   coefficient vectors to represent multi-precision natural numbers.

-- SHIFTS AND NORMALIZATIONS :

  procedure Digits_to_Left ( n : in out Array_of_Naturals; k : in natural32 );

  -- DESCRIPTION :
  --   Shifts all digits of n to the left by k places,
  --   corresponding to a multiplication by 10^k.

  procedure Shift_Left ( n : in out Array_of_Naturals; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts n until the leading coefficient is larger than Base/Radix.
  --   The number ns on return equals the number of shifts.

  procedure Shift_Right ( n : in out Array_of_Naturals; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts n until the last digit in the least significant coefficient
  --   of n is different from zero.
  --   The number ns on return equals the number of shifts.

-- OPERATIONS :

  procedure Add ( n1 : in out Array_of_Naturals; n2 : in natural64 );

  -- DESCRIPTION :
  --   Adds the number n2 to the coefficient vector n1.

  -- REQUIRED : n1 is long enough to hold the result.
  --   No exception is generated when the carry over is nonzero and
  --   n1(n1'last) can no longer contain it.

  procedure Add ( n1 : in out Array_of_Naturals;
                  n2 : in Array_of_Naturals );

  -- DESCRIPTION :
  --   Adds the coefficient vector representation n2 to n1.

  -- REQUIRED : n1 > n2 and n1 is long enough to hold the result.

  procedure Sub ( n1 : in out Array_of_Naturals; n2 : in Array_of_Naturals );

  -- DESCRIPTION :
  --   Stores the coefficient vector representation of n1-n2 in n1.

  -- REQUIRED : n1 >= n2.  No exception is generated if this is not the case.

  procedure Mul_Fact ( n : in out Array_of_Naturals; f : in natural64 );

  -- DESCRIPTION :
  --   Multiplies the number n with a factor f that is smaller than the radix.
  --   The result is stored in n.
  
  -- REQUIRED : n is large enough to hold the result.

  procedure Acc_Add ( n : in out Array_of_Naturals; m1,m0 : in natural64;
                      k : in natural32; c : in natural64 );

  -- DESCRIPTION :
  --   Accumulated addition: m0 is added to n at position k, and m1 is
  --   added, eventually with a carry over to n at position k+1.
  --   The eventual carry over of the that addition goes to position k+2.
  --   This routine is auxiliary to the multiplication routines.

  procedure Mul ( n1 : in out Array_of_Naturals; n2 : in natural64 );

  -- DESCRIPTION 
  --   Stores the result of n1*n2 in n1.

  -- REQUIRED : n1 /= 0 /= n2.
  --   n1 contains zeroes in its two last entries for carry over.

  function Mul ( n1,n2 : Array_of_Naturals ) return Array_of_Naturals;

  -- DESCRIPTION :
  --   Returns the results of the product n1*n2.

  -- REQUIRED : n1'last >= n2'last.

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural64 );

  -- DESCRIPTION :
  --   Computes the quotient in the coefficient vector representation
  --   of the division n1/n2.

  -- REQUIRED : n2 /= 0.

  procedure Small_Div ( n1 : in Array_of_Naturals; n2 : in natural64;
                        q : out Array_of_Naturals; r : out natural64 );

  -- DESCRIPTION :
  --   n1 = n2*q + r, only applies when n2 <= radix.

  -- REQUIRED : q'range = n1'range and n2 /= 0 /= n1 and n1 > n2.

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural64;
                        r : out natural64 );

  -- DESCRIPTION :
  --   The quotient of the division n1/n2 is stored in n1.
  --   Same restrictions as with the other Small_Div apply.

  procedure Big_Div ( n1 : in Array_of_Naturals; n2 : in natural64;
                      q : out Array_of_Naturals; r : out natural64 );

  -- DESCRIPTION :
  --   This is the division for n2 > radix on the coefficient vector
  --   representation of a multi-precision natural number n1.

  -- REQUIRED : q'range = n1'range; n2 /= 0 /= n1 and n1 > n2.

  procedure Big_Div ( n1 : in out Array_of_Naturals; n2 : in natural64;
                      r : out natural64 );

  -- DESCRIPTION :
  --   The quotient of the division n1/n2 is stored in n1.
  --   Same restrictions as with the other Big_Div apply.

  procedure Div ( n1,n2 : in Array_of_Naturals;
                  q,r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   Performs the Div operation on the coefficient vector representations
  --   of multi-precision natural numbers.

  -- REQUIRED :
  --   q'range contains n1'range, r'range contains n2'range.
  --   Trivial cases have been taken care of, in particular:
  --   n1 > n2 and n2'last > 0 with n2(n2'last) /= 0.

  procedure Div ( n1 : in out Array_of_Naturals; n2 : in Array_of_Naturals;
                  r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   The resulting quotient of the division n1/n2 is stored in n1.
  --   Same restrictions as the above Div operation apply.

end Multprec_Natural64_Coefficients;
