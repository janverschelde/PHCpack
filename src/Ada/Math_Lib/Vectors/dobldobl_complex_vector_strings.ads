with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;

package DoblDobl_Complex_Vector_Strings is

-- DESCRIPTION :
--   Exports functions to write a double double complex vector to string
--   and to parse a vector from a string.
--   The simple format of a complex number is followed:
--   a complex number is a sequence of two double doubles,
--   separated by two spaces.

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

end DoblDobl_Complex_Vector_Strings; 
