with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with Multprec_Complex_numbers;

package Strings_and_Numbers is

-- DESCRIPTION :
--   The functions in this package write numbers to strings
--   and parse strings back into numbers.

  function Is_Imag ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean;
  function Is_Imag ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if the number is pure imaginary, false otherwise.

  function Is_Real ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean;
  function Is_Real ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if the complex number is real, false otherwise.

  function Is_Integer ( f : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the floating number is an integer, false otherwise.

  function Is_Unit ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean;
  function Is_Unit ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if c is +1 or -1.

  function Sign ( c : Standard_Complex_Numbers.Complex_Number )
                return integer;
  function Sign ( c : Multprec_Complex_Numbers.Complex_Number )
                return integer;

  -- DESCRIPTION :
  --   This routine only make sense when Is_Unit(c).
  --   It returns +1 or -1, depending whether c = +1 or -1.

  function Truncate ( f : double_float ) return integer32;

  -- DESCRIPTION :
  --   Truncates the double float into a 32-bit integer.

  function Trim_Zeros ( s : string ) return string;

  -- DESCRIPTION :
  --   Returns the string s with trailing zeros removed.

  function Convert ( f : double_float ) return string;
  function Convert ( f : Floating_Number ) return string;

  -- DESCRIPTION :
  --   Converts a positive float into fixed point format string.

  function Convert
             ( c : Standard_Complex_Numbers.Complex_Number ) return string;
  function Convert
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string;

  -- DESCRIPTION :
  --   Returns a string with a complex number in the format
  --   "( re + im*i)", used for writing coefficients of polynomials.

  function Signed_Constant
              ( c : Standard_Complex_Numbers.Complex_Number ) return string;
  function Signed_Constant
              ( c : Multprec_Complex_Numbers.Complex_Number ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the constant term of
  --   a polynomial, with its sign.

  function Unsigned_Constant
              ( c : Standard_Complex_Numbers.Complex_Number ) return string;
  function Unsigned_Constant
              ( c : Multprec_Complex_Numbers.Complex_Number ) return string;

  -- DESCRIPTION :
  --   This routine converts the constant term of the polynomial into
  --   a string.  A positive integer is returned without the sign.

  function Signed_Coefficient
              ( c : Standard_Complex_Numbers.Complex_Number ) return string;
  function Signed_Coefficient
              ( c : Multprec_Complex_Numbers.Complex_Number ) return string;

  -- DESCRIPTION :
  --   Returns the coefficient with its sign.  
  --   The +1 and -1 are treated differently.

  function Unsigned_Coefficient
              ( c : Standard_Complex_Numbers.Complex_Number ) return string;
  function Unsigned_Coefficient
              ( c : Multprec_Complex_Numbers.Complex_Number ) return string;

  -- DESCRIPTION :
  --   Returns the coefficient without sign, except if it is negative.

end Strings_and_Numbers;
