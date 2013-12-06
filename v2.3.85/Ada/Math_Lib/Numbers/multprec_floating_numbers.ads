with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;

package Multprec_Floating_Numbers is

-- DESCRIPTION :
--   This package provides operations to deal with floating-point numbers
--   of arbitrary length.  The size of the fraction remains fixed during
--   computations, but can be adjusted.  The exponent has no fixed size.
--   Mixed-precision arithmetic is fully supported: the size of the result
--   in any binary operation is the maximal size of the operands.

  type Floating_Number is private;

-- CONSTRUCTORS :

  function Create ( i : integer ) return Floating_Number;
  function Create ( n : natural32 ) return Floating_Number;
  function Create ( i : integer32 ) return Floating_Number;
  function Create ( i : Integer_Number ) return Floating_Number;
  function Create ( f : double_float ) return Floating_Number;

  function Create ( fraction,exponent : integer32 ) return Floating_Number;
  function Create ( fraction : Integer_Number;
                    exponent : integer32 ) return Floating_Number;
  function Create ( fraction : integer32;
                    exponent : Integer_Number ) return Floating_Number;
  function Create ( fraction,exponent : Integer_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   A floating point number consists of a fraction and exponent.

-- SELECTORS :

  function Fraction ( f : Floating_Number ) return Integer_Number;
  function Exponent ( f : Floating_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the fraction and the exponent.  Warning: sharing!
  --   We have that f = Fraction(f)*10**Exponent(f).

  function Size_Fraction ( f : Floating_Number ) return natural32;
  function Size_Exponent ( f : Floating_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the length of the coefficient vector of the representations
  --   of fraction and exponent as integer numbers.
  --   Size_Fraction(f) = 1 corresponds to the precision of the standard
  --   double floating-point numbers.

  function Decimal_Places_Fraction ( f : Floating_Number ) return natural32;
  function Decimal_Places_Exponent ( f : Floating_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places in fraction and exponent.

  function Decimal_to_Size ( deci : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the smallest size to obtain the given number of decimal places.
  --   This is useful for the "Set_Size" operation.

  function Size_to_Decimal ( size : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places corresponding to the size.

  function AbsVal ( f : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the absolute value of f.  There is no sharing with f.

  function Positive ( f : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the fraction of f is positive.

  function Negative ( f : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the fraction of f is negative.

  function Sign ( f : Floating_Number ) return integer32;

  -- DESCRIPTION :
  --   Returns the sign of the fraction of f.

  function Truncate_to_Nearest_Integer
             ( f : Floating_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the integer nearest to f after truncation.

  -- REQUIRED : Exponent(f) fits into a 32-bit integer.

-- COMPARISON AND COPYING :

  function Equal ( f1 : Floating_Number; f2 : double_float ) return boolean;
  function Equal ( f1,f2 : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if f1 and f2 are equal, false otherwise.

  function "<" ( f1 : double_float; f2 : Floating_Number ) return boolean;
  function "<" ( f1 : Floating_Number; f2 : double_float ) return boolean;
  function "<" ( f1,f2 : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if f1 < f2, false otherwise.

  function ">" ( f1 : double_float; f2 : Floating_Number ) return boolean;
  function ">" ( f1 : Floating_Number; f2 : double_float ) return boolean;
  function ">" ( f1,f2 : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if f1 > f2, false otherwise;

  procedure Copy ( f1 : in Floating_Number; f2 : in out Floating_Number );

  -- DESCRIPTION :
  --   Clears f2 and copies the content of f1 into f2.
  --   Note that f2 := f1 leads to data sharing.

-- EXPANDING and SHORTENING THE MANTISSA :

  function  Expand ( f : Floating_Number; k : natural32 )
                   return Floating_Number;
  procedure Expand ( f : in out Floating_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Expands the fraction by adding k coefficients.

  function  Round ( f : Floating_Number ) return double_float;
  function  Round ( f : Floating_Number; k : natural32 ) return Floating_Number;
  procedure Round ( f : in out Floating_Number; k : in natural32 );

  function  Trunc ( f : Floating_Number ) return double_float;
  function  Trunc ( f : Floating_Number; k : natural32 ) return Floating_Number;
  procedure Trunc ( f : in out Floating_Number; k : in natural32 );

  -- DESCIRPTION :
  --   Shortens the fraction by removing k coefficients, either by rouding
  --   or by truncation of the least significant coefficients.

  procedure Set_Size ( f : in out Floating_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Sets the size of the floating number f to k, either by rounding
  --   or expanding its current size.

-- ARITHMETIC OPERATIONS as functions (no data sharing) :
--   The size of the result is the maximal size of the operands.

  function "+" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number;
  function "+" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number;
  function "+" ( f1,f2 : Floating_Number ) return Floating_Number;

  function "+" ( f : Floating_Number ) return Floating_Number;  -- copy of f
  function "-" ( f : Floating_Number ) return Floating_Number;

  function "-" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number;
  function "-" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number;
  function "-" ( f1,f2 : Floating_Number ) return Floating_Number;

  function "*" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number;
  function "*" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number;
  function "*" ( f1,f2 : Floating_Number ) return Floating_Number;

  function "**" ( f : double_float; n : Natural_Number ) return Floating_Number;
  function "**" ( f : Floating_Number; n : Natural_Number )
                return Floating_Number;

  function "**" ( f : Floating_Number; i : integer32 ) return Floating_Number;
  function "**" ( f : double_float; i : Integer_Number ) return Floating_Number;
  function "**" ( f : Floating_Number; i : Integer_Number )
                return Floating_Number;

  function "/" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number;
  function "/" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number;
  function "/" ( f1,f2 : Floating_Number ) return Floating_Number;

-- ARITHMETIC OPERATIONS as procedures for memory management :

  procedure Add ( f1 : in out Floating_Number; f2 : in double_float );  -- "+"
  procedure Add ( f1 : in out Floating_Number; f2 : in Floating_Number );

  procedure Sub ( f1 : in out Floating_Number; f2 : in double_float );  -- "-"
  procedure Sub ( f1 : in out Floating_Number; f2 : in Floating_Number );
  procedure Min ( f : in out Floating_Number );

  procedure Mul ( f1 : in out Floating_Number; f2 : in double_float );  -- "*"
  procedure Mul ( f1 : in out Floating_Number; f2 : in Floating_Number );

  procedure Div ( f1 : in out Floating_Number; f2 : in double_float );  -- "/"
  procedure Div ( f1 : in out Floating_Number; f2 : in Floating_Number );

-- DESTRUCTOR :

  procedure Clear ( f : in out Floating_Number );

  -- DESCRIPTION :
  --   Deallocations of the occupied memory for f.

private

  type Floating_Number is record
    fraction,exponent : Integer_Number;
  end record;

end Multprec_Floating_Numbers;
