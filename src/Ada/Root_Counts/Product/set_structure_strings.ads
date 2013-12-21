with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;

package Set_Structure_Strings is

-- DESCRIPTION :
--   Writing and parsing of set structures to and from strings.

-- WRITE STRING REPRESENTATIONS :

  function to_String ( i,j : in natural32 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the j-th set
  --   of the set structure for the j-th polynomial.
  --   The format is '{' followed by the strings for the symbols,
  --   separated by on space, terminated by '}'.

  function to_String ( i : in natural32 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the set structure
  --   for the i-th polynomial, which consists of the string
  --   representations for all sets in the structure for the i-th polynomial.

  function to_String return string;

  -- DESCRIPTION :
  --   Returns the string representation of the entire set structure.
  --   This string consists in the string representations of the
  --   set structures for the polynomials, separated by semicolons ';'.

-- PARSE STRING REPRESENTATIONS :

  function Number_of_Semicolons ( s : string ) return natural32;

  -- DESCRIPTION :
  --   Counts the number of semicolons in the string.

  function Number_of_Sets
             ( s : string ) return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..Number_of_Semicolons(s)
  --   with the i-th entry in the vector on return equal to the
  --   number of sets for the set structure of the i-th polynomial.

  procedure Parse ( s : in string );

  -- DESCRIPTION :
  --   Parses the string s for a set structure, splitting the string
  --   on the semicolons and filling the set structure package data
  --   with the set structures of the polynomials.

  procedure Parse ( s : in string; i : in natural32 );

  -- DESCRIPTION :
  --   Parses the string s into the set structure for the i-th polynomial.

end Set_Structure_Strings;
