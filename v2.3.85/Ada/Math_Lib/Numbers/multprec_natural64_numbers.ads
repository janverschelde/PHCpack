with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural64_Coefficients;    use Multprec_Natural64_Coefficients;

package Multprec_Natural64_Numbers is

-- DESCRIPTION :
--   This package allows to manipulate natural numbers of arbitrary length,
--   based on 64-bit integer arithmetic.

-- DATA STRUCTURE :

  type Natural_Number is private;

-- CREATORS :

  function Create ( n : natural64 ) return Array_of_Naturals;

  -- DESCRIPTION :
  --   Returns the representation of the natural number as coefficient vector.

  function Create ( n : natural64 ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the representation of n as a natural number.

  function Create ( n : Array_of_Naturals ) return Natural_Number;

  -- DESCRIPTION :
  --   Creates a number from the coefficients in n.

  -- REQUIRED : n'range = 0..n'last and n(i) < Base, for i in n'range.

  function Create ( n : Natural_Number ) return natural64;

  -- DESCRIPTION :
  --   Returns the representation of n as a standard natural.

  -- REQUIRED : n < Base.

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

  function Empty ( n : Natural_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the number has not been created yet, or when it has
  --   been destroyed by the operation Clear; otherwise false is returned.
  --   An empty number is considered as zero.

  function Size ( n : Natural_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the last entry in the coefficient representation.

  function Decimal_Places ( n : natural64 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places n occupies, Decimal_Places(0) = 0.

  function Decimal_Places ( n : Natural_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places n occupies, Decimal_Places(0) = 0.
  --   Since n can be arbitrarily large, also the number of return should
  --   have no constraints on its size.  However, the current implementation
  --   would not support the manipulation of numbers whose size causes this
  --   function to crash.

  function Coefficient ( n : Natural_Number; i : natural32 ) return natural64;

  -- DESCRIPTION :
  --   Returns the ith entry in the coefficient representation.

  function Coefficients ( n : Natural_Number ) return Array_of_Naturals;

  -- DESCRIPTION :
  --   Returns the coefficient representation of n, of range 0..Size(n).
  --   The number n equals then the sum of Coefficient(n,i)*Base**i,
  --   for i in 0..Size(n).

-- COMPARISON AND COPYING :

  function Equal ( n1 : Natural_Number; n2 : natural64 ) return boolean;
  function Equal ( n1,n2 : Natural_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true when both numbers n1 and n2 are equal, false otherwise.

  function "<" ( n1 : Natural_Number; n2 : natural64 ) return boolean;
  function "<" ( n1 : natural64; n2 : Natural_Number ) return boolean;
  function "<" ( n1,n2 : Natural_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if n1 < n2, false otherwise.

  function ">" ( n1 : Natural_Number; n2 : natural64 ) return boolean;
  function ">" ( n1 : natural64; n2 : Natural_Number ) return boolean;
  function ">" ( n1,n2 : Natural_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if n1 > n2, false otherwise.

  procedure Copy ( n1 : in natural64; n2 : in out Natural_Number );
  procedure Copy ( n1 : in Natural_Number; n2 : in out Natural_Number );

  -- DESCRIPTION :
  --   Clears n2 and makes a copy of n1 to be equal to n2.
  --   Note that n2 := n1 leads to data sharing.

-- SHIFTS :

  procedure Shift_Left ( n : in out Natural_Number; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts n until the leading coefficient is larger than Base/Radix.
  --   The number ns on return equals the number of shifts.

  procedure Shift_Right ( n : in out Natural_Number; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts n until the last digit in the least significant coefficient
  --   of n is different from zero.
  --   The number ns on return equals the number of shifts.

-- ARITHMETIC OPERATIONS as functions (no data sharing) :
--   Note that n1 >= n2 is required for subtraction, and n2 /= 0 for division.
--   The unary "-" operations have been added to make it ring-like.

  function "+" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number;
  function "+" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number;
  function "+" ( n1,n2 : Natural_Number ) return Natural_Number;

  function "-" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number;
  function "-" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number;
  function "-" ( n1,n2 : Natural_Number ) return Natural_Number;

  function "+" ( n : Natural_Number ) return Natural_Number;   -- copies n
  function "-" ( n : Natural_Number ) return Natural_Number;   -- copies n

  function "*" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number;
  function "*" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number;
  function "*" ( n1,n2 : Natural_Number ) return Natural_Number;

  function "**" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number;
  function "**" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number;
  function "**" ( n1,n2 : Natural_Number ) return Natural_Number;

  function "/" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number;
  function "/" ( n1 : natural64; n2 : Natural_Number ) return natural64;
  function "/" ( n1,n2 : Natural_Number ) return Natural_Number;

  function Rmd ( n1 : Natural_Number; n2 : natural64 ) return natural64;
  function Rmd ( n1 : natural64; n2 : Natural_Number ) return natural64;
  function Rmd ( n1,n2 : Natural_Number ) return Natural_Number;

-- SPECIAL ARITHMETIC OPERATION : times radix = shift

  procedure Mul_Radix ( n : in out Natural_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Multiplies n with the radix k times, which corresponds to a shift
  --   of n with k places to the left.

-- ARITHMETIC OPERATIONS as procedures for memory management :

  procedure Add ( n1 : in out Natural_Number; n2 : in natural64 );    -- "+"
  procedure Add ( n1 : in out Natural_Number; n2 : in Natural_Number );

  procedure Sub ( n1 : in out Natural_Number; n2 : in natural64 );    -- "-"
  procedure Sub ( n1 : in out Natural_Number; n2 : in Natural_Number );
  procedure Min ( n : in out Natural_Number );

  procedure Mul ( n1 : in out Natural_Number; n2 : in natural64 );    -- "*"
  procedure Mul ( n1 : in out Natural_Number; n2 : in Natural_Number );

  procedure Div ( n1 : in out Natural_Number; n2 : in natural64 );    -- "/"
  procedure Div ( n1 : in out Natural_Number; n2 : in Natural_Number );

  procedure Div ( n1 : in Natural_Number; n2 : in natural64;  -- n1 = n2*q+r
                  q : out Natural_Number; r : out natural64 );
  procedure Div ( n1 : in out Natural_Number; n2 : in natural64;
                  r : out natural64 );
  procedure Div ( n1,n2 : in Natural_Number; q,r : out Natural_Number );
  procedure Div ( n1 : in out Natural_Number; n2 : in Natural_Number;
                  r : out Natural_Number );

-- DESTRUCTOR :

  procedure Clear ( n : in out Natural_Number );

  -- DESCRIPTION :
  --   Deallocation of the memory space.  Empty(n) is true on return.

private
 
  type Natural_Number_Rep;
  type Natural_Number is access Natural_Number_Rep;

end Multprec_Natural64_Numbers;
