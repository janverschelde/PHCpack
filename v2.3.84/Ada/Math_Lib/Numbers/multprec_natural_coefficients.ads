with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Multprec_Natural_Coefficients is

-- DESCRIPTION :
--   This package provides manipulations on coefficient vectors of natural
--   numbers to support multi-precision arithmetic.
--   The coefficient vector c represents a multi-precision natural number
--   as the sum of c(i)*Base^i, with c(i) < Base, for i in 0..c'last.
--   The operations in this package are provided for efficiency.
--   No automatic tests are built in to check the requirements.
--   Routines that call these operations must make sure the requirements
--   for the operations are met.

-- DATASTRUCTURE :

  type Array_of_Naturals is array ( natural32 range <> ) of natural32;

-- SELECTORS :

  function Radix return natural32;

  -- DESCRIPTION :
  --   Returns the radix of the number representation.

  function Base return natural32;

  -- DESCRIPTION :
  --   Returns the base = Radix^Exponent of the number representation.

  function Exponent return natural32;

  -- DESCRIPTION :
  --   Returns the exponent which determines the size of the base.

  function Size_of_Coefficient ( n : natural32 ) return natural32;

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

  procedure Add ( n1 : in out Array_of_Naturals; n2 : in natural32 );

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

  procedure Mul_Fact ( n : in out Array_of_Naturals; f : in natural32 );

  -- DESCRIPTION :
  --   Multiplies the number n with a factor f that is smaller than the radix.
  --   The result is stored in n.
  --   Note that the opposite operation is achieved by Small_Div below.
  
  -- REQUIRED : n is large enough to hold the result.

  procedure Short_Mul ( n1 : in out Array_of_Naturals; n2 : in natural32 );

  -- DESCRIPTION 
  --   Stores the result of n1*n2 in n1.

  -- REQUIRED : n1 /= 0 /= n2, n2 < Base,
  --   n1 contains zeroes in its two last entries for carry over.

  function Mul ( n1,n2 : Array_of_Naturals ) return Array_of_Naturals;

  -- DESCRIPTION :
  --   Returns the results of the product n1*n2.

  -- REQUIRED : n1'last >= n2'last.

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural32 );

  -- DESCRIPTION :
  --   Computes the quotient in the coefficient vector representation
  --   of the division n1/n2.

  -- REQUIRED : 0 < n2 <= Radix.

  procedure Small_Div ( n1 : in Array_of_Naturals; n2 : in natural32;
                        q : out Array_of_Naturals; r : out natural32 );

  -- DESCRIPTION :
  --   n1 = n2*q + r, only applies when n2 <= Radix.

  -- REQUIRED :
  --   q'range = n1'range and n2 /= 0 /= n1 and n1 > n2, n2 <= Radix.

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural32;
                        r : out natural32 );

  -- DESCRIPTION :
  --   The quotient of the division n1/n2 is stored in n1.
  --   Same restrictions as with the other Small_Div apply.

  procedure Short_Div ( n1 : in Array_of_Naturals; n2 : in natural32;
                        q : out Array_of_Naturals; r : out natural32 );

  -- DESCRIPTION :
  --   The division runs much faster when the divisor is a coefficient.

  -- REQUIRED : 0 < n2 < Base and q'range = n1'range.

  procedure Short_Div ( n1 : in out Array_of_Naturals; n2 : in natural32;
                        r : out natural32 );

  -- DESCRIPTION :
  --   The quotient is stored in n1.

  -- REQUIRED : 0 < n2 < Base.

  procedure Short_Div2 ( n1,n2 : in Array_of_Naturals;
                         q,r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   n1 = n2*q + r in case n2 = n2(1)*Base + n2(0).

  -- REQUIRED :
  --   n2 > 0 and n2'last = 1 or n2(i) = 0 for i > 1,
  --   for ranges: q'last >= n1'last and r'last >= n2'last.

  procedure Short_Div2 ( n1 : in out Array_of_Naturals;
                         n2 : in Array_of_Naturals;
                         r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   Overwrites n1 with the quotient q, n1 = q*n2 + r.
  --   The same constraints apply as in the other Short_Div2.

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

  procedure Div2 ( n1,n2 : in Array_of_Naturals;
                   q,r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   Returns in q and r the quotient and remainder of division of n1 by n2,
  --   i.e.: n1 = n2*q + r on return.

  -- REQUIRED : n2'last > 1 and n2(n2'last) /= 0.
  --   Otherwise the Short_Div's should be applied.

  procedure Div2 ( n1 : in out Array_of_Naturals; n2 : in Array_of_Naturals;
                   r : out Array_of_Naturals );

  -- DESCRIPTION :
  --   Overwrites the n1 with the quotient, so on return n1 = n1*n2 + r.
  --   The same constraints as with the other Div2 hold.

end Multprec_Natural_Coefficients;
