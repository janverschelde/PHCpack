with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;

package Sets_of_Unknowns_Strings is

-- DESCRIPTION :
--   This package provides operations to write sets of unknowns to strings
--   and to parse strings into sets of unknowns.

  function to_String ( element : natural32 ) return string;

  -- DESCRIPTION :
  --   Writes the element of a set to string.

  function to_String ( s : Set ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the set.

  function Parse ( s : string; n : natural32 ) return Set;

  -- DESCRIPTION :
  --   Parses the string into a set of n unknowns.

  -- REQUIRED :
  --   The string must start with { and end with }.
  --   Variables must be separated by spaces.
  --   If the symbol table is empty, then the first character
  --   of each variable name must be 'x' and be followed by the number,
  --   for example 'x18' to denote the 18-th variable.
  --   If the symbol table is not empty, then the variable names in 
  --   the string s must be present in the symbol table.

end Sets_of_Unknowns_Strings;
