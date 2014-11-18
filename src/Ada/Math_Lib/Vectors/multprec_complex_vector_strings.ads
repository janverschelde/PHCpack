with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;

package Multprec_Complex_Vector_Strings is

-- DESCRIPTION :
--   Exports functions to write a multiprecision complex vector 
--   to a string and to parse a string into a multprecision complex vector.
--   The simple format of a complex number is followed:
--   a complex number is a sequence of two multiprecision floating-point
--   numbers (real and imaginary part), separated by two spaces.

  function Write ( v : Vector ) return string;

  -- DESCRIPTION :
  --   Writes the numbers in the vector v to string in the full
  --   scientific format, separated by newline symbols.
  --   Real and imaginary parts of the complex numbers are
  --   separated by spaces.

  function Parse ( s : string ) return Vector;

  -- DESCRIPTION :
  --   Parses the numbers separated by newline symbols in the string
  --   to the vector on return.

end Multprec_Complex_Vector_Strings; 
