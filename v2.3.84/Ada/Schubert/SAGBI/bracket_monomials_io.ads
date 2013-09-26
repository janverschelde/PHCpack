with text_io;                          use text_io;
with Brackets,Bracket_Monomials;       use Brackets,Bracket_Monomials;

package Bracket_Monomials_io is

-- DESCRIPTION :
--   This package provides input/output operations for bracket monomials.

  procedure get ( b : out Link_to_Bracket );
  procedure get ( file : in file_type; b : out Link_to_Bracket );

  -- DESCRIPTION :
  --   Reads a natural vector, enclosed by square brackets.
  --   The numbers in the bracket need to be separated by spaces,
  --   spaces before or after the brackets are ignored.

  procedure get ( bm : out Bracket_Monomial );
  procedure get ( file : in file_type;bm : out Bracket_Monomial );

  -- DESCRIPTION :
  --   Reads a bracket monomial in the format as in [1 2 3]^3*[2 3 4];

  procedure put ( bm : in Bracket_Monomial );
  procedure put ( file : in file_type; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Writes a bracket monomial on standard output or on file if
  --   this file is specified as argument.

end Bracket_Monomials_io;
