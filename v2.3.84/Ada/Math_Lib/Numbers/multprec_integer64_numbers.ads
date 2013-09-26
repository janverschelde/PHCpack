with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Natural64_Coefficients;    use Multprec_Natural64_Coefficients;
with Multprec_Natural64_Numbers;         use Multprec_Natural64_Numbers;

package Multprec_Integer64_Numbers is

-- DESCRIPTION :
--   This package allows to manipulate integer numbers of arbitrary length.
--   The operations are such that this is a faithful representation of the
--   integers as an Euclidean domain.

-- DATA STRUCTURES :

  type Integer_Number is private;

-- CREATORS :

  function Shallow_Create ( n : Natural_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   This create returns an integer number that shares the n.
  --   If n changes, also the result of this function will change.
  --   To clear the result of this function, use Shallow_Clear.

  function Create ( n : Natural_Number ) return Integer_Number;
  function Create ( n : Array_of_Naturals ) return Integer_Number;
  function Create ( i : integer ) return Integer_Number;
  function Create64 ( i : integer64 ) return Integer_Number;
  function Create32 ( i : integer32 ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the representation of n or i as a integer number.

  function Convert ( n : Natural_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   This operation has the same effect as create, except that
  --   n and i := Convert(n), share the same number: if n changes
  --   after conversion, i will change as well, and vice versa.
  --   This operation can be useful for integer/natural arithmetic.

  function Create ( i : Integer_Number ) return integer64;

  -- DESCRIPTION :
  --   Returns the representation of i as a standard integer.

  -- REQUIRED : i < Natural_Numbers.Basis.

-- SELECTORS :

  function Empty ( i : Integer_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the number has not been created yet, or when it has
  --   been destroyed by the operation Clear; otherwise false is returned.
  --   An empty number is considered as zero.

  function Size ( i : Integer_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the last entry in the coefficient representation.

  function Coefficient ( i : Integer_Number; k : natural32 ) return natural64;

  -- DESCRIPTION :
  --   Returns the kth component in the coefficient representation of i.

  function Coefficients ( i : Integer_Number ) return Array_of_Naturals;

  -- DESCRIPTION :
  --   Returns the coefficient representation of the unsigned integer.

  function Decimal_Places ( i : Integer_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places of i.
  --   The number of decimal places of 0 is 0.

  function Positive ( i : Integer_Number ) return boolean;
  function Negative ( i : Integer_Number ) return boolean;

  -- DESCRIPTION :
  --   If i > 0, then Positive(i) and not Negative(i).
  --   If i < 0, then not Positive(i) and Negative(i).
  --   For i = 0, sign switching by Min(i) can give either +0 or -0.

  function Sign ( i : Integer_Number ) return integer32;

  -- DESCRIPTION :
  --   Returns +1,-1, or 0, depending whether i > 0, i < 0, or i = 0.

  function Unsigned ( i : Integer_Number ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the unsigned integer number.
  --   Note that this is a pointer with data sharing!

-- COMPARISON AND COPYING :

  function Equal ( i1 : Integer_Number; i2 : integer64 ) return boolean;
  function Equal ( i1,i2 : Integer_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true when both numbers i1 and i2 are equal, false otherwise.

  function "<" ( i1 : Integer_Number; i2 : integer64 ) return boolean;
  function "<" ( i1 : integer64; i2 : Integer_Number ) return boolean;
  function "<" ( i1,i2 : Integer_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if i1 < i2, false otherwise.

  function ">" ( i1 : Integer_Number; i2 : integer64 ) return boolean;
  function ">" ( i1 : integer64; i2 : Integer_Number ) return boolean;
  function ">" ( i1,i2 : Integer_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if i1 > i2, false otherwise.

  procedure Copy ( i1 : in integer64; i2 : in out Integer_Number );
  procedure Copy ( i1 : in Integer_Number; i2 : in out Integer_Number );

  -- DESCRIPTION :
  --   Clears i2 and makes a copy of i1 to be equal to i2.
  --   Note that i2 := i1 leads to data sharing.

-- SHIFTS :

  procedure Shift_Left ( i : in out Integer_Number; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts i until the leading coefficient is larger than Base/Radix.
  --   The number ns on return equals the number of shifts.

  procedure Shift_Right ( i : in out Integer_Number; ns : out natural32 );

  -- DESCRIPTION :
  --   Shifts i until the last digit in the least significant coefficient
  --   of i is different from zero.
  --   The number ns on return equals the number of shifts.

-- ARITHMETIC OPERATIONS as functions (no data sharing) :

  function "+" ( i1 : Integer_Number; i2 : integer64 ) return Integer_Number;
  function "+" ( i1 : integer64; i2 : Integer_Number ) return Integer_Number;
  function "+" ( i1,i2 : Integer_Number ) return Integer_Number;

  function "+" ( i : Integer_Number ) return Integer_Number;  -- copies i
  function "-" ( i : Integer_Number ) return Integer_Number;

  function "-" ( i1 : Integer_Number; i2 : integer64 ) return Integer_Number;
  function "-" ( i1 : integer64; i2 : Integer_Number ) return Integer_Number;
  function "-" ( i1,i2 : Integer_Number ) return Integer_Number;

  function "*" ( i1 : Integer_Number; i2 : integer64 ) return Integer_Number;
  function "*" ( i1 : integer64; i2 : Integer_Number ) return Integer_Number;
  function "*" ( i1,i2 : Integer_Number ) return Integer_Number;

  function "**" ( i : Integer_Number; n : natural64 ) return Integer_Number;
  function "**" ( i : integer64; n : Natural_Number ) return Integer_Number;
  function "**" ( i : Integer_Number; n : Natural_Number )
                return Integer_Number;

  function "/" ( i1 : Integer_Number; i2 : integer64 ) return Integer_Number;
  function "/" ( i1 : integer64; i2 : Integer_Number ) return integer64;
  function "/" ( i1,i2 : Integer_Number ) return Integer_Number;

  function Rmd ( i1 : Integer_Number; i2 : integer64 ) return integer64;
  function Rmd ( i1 : integer64; i2 : Integer_Number ) return integer64;
  function Rmd ( i1,i2 : Integer_Number ) return Integer_Number;

-- SPECIAL ARITHMETIC OPERATION : multiply with radix = shift left

  procedure Mul_Radix ( i : in out Integer_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Multiplies i with the radix k times.

-- ARITHMETIC OPERATIONS as procedures for memory management :

  procedure Add ( i1 : in out Integer_Number; i2 : in integer64 );    -- "+"
  procedure Add ( i1 : in out Integer_Number; i2 : in Integer_Number );

  procedure Min ( i : in out Integer_Number );
  procedure Set_Min ( i : in out Integer_Number ); -- allows also -0

  procedure Sub ( i1 : in out Integer_Number; i2 : in integer64 );    -- "-"
  procedure Sub ( i1 : in out Integer_Number; i2 : in Integer_Number );

  procedure Mul ( i1 : in out Integer_Number; i2 : in integer64 );    -- "*"
  procedure Mul ( i1 : in out Integer_Number; i2 : in Integer_Number );

  procedure Rmd ( i1 : in out Integer_Number; i2 : in integer64 );
  procedure Rmd ( i1 : in out Integer_Number; i2 : in Integer_Number );

  procedure Div ( i1 : in out Integer_Number; i2 : in integer64 );    -- "/"
  procedure Div ( i1 : in out Integer_Number; i2 : in Integer_Number );

  procedure Div ( i1 : in Integer_Number; i2 : in integer64;  -- i1 = i2*q+r
                  q : out Integer_Number; r : out integer64 );
  procedure Div ( i1 : in out Integer_Number; i2 : in integer64;
                  r : out integer64 );
  procedure Div ( i1,i2 : in Integer_Number; q,r : out Integer_Number );
  procedure Div ( i1 : in out Integer_Number; i2 : in Integer_Number;
                  r : out Integer_Number );

-- DESTRUCTOR :

  procedure Shallow_Clear ( i : in out Integer_Number );

  -- DESCRIPTION :
  --   Clears only the memory needed to store i, not the Natural_Number
  --   that i contains.

  procedure Clear ( i : in out Integer_Number );

  -- DESCRIPTION :
  --   Deallocation of the memory space.  Empty(i) is true on return.

private
 
  type Integer_Number_Rep;
  type Integer_Number is access Integer_Number_Rep;

end Multprec_Integer64_Numbers;
