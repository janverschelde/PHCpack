with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;

package Standard_Complex_Vector_Strings is

-- DESCRIPTION :
--   Exports functions to write a standard complex vector to string
--   and to parse a vector from a string.
--   The simple format of a complex number is followed:
--   a complex number is a sequence of two double floats,
--   separated by two spaces.

  function Write ( v : Vector ) return string;

  -- DESCRIPTION :
  --   Writes the numbers in the vector v to string in the full
  --   scientific format, separated by newline symbols.
  --   Real and imaginary parts of the complex numbers are
  --   separated by spaces.

  function Count_Linefeeds ( s : string ) return integer32;

  -- DESCRIPTION :
  --   Counts the number of line feeds in the string.

  function Next_Linefeed ( s : string ) return integer;

  -- DESCRIPTION :
  --   Returns the position in the string where the next linefeed occurs.
  --   If there is no linefeed in s, then the number on return is
  --   larger than s'last.

  function Parse ( s : string ) return Vector;

  -- DESCRIPTION :
  --   Parses the numbers separated by newline symbols in the string
  --   to the vector on return.

end Standard_Complex_Vector_Strings; 
