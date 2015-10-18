with text_io;                        use text_io;
with QuadDobl_Bracket_Polynomials;   use QuadDobl_Bracket_Polynomials;

package QuadDobl_Bracket_Polynomials_io is

-- DESCRIPTION :
--   This package provides output operations for bracket polynomials,
--   with complex coefficients in double double precision.

  procedure put ( t : in Bracket_Term );
  procedure put ( file : in file_type; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Writes the term on standard output or on file,
  --   depending whether a file is an argument of the procedure.

  procedure put ( p : in Bracket_Polynomial );
  procedure put ( file : in file_type; p : in Bracket_Polynomial );

  -- DESCRIPTION :
  --   Writes the polynomial on standard output or on file,
  --   depending whether a file is an argument of the procedure.

  procedure put_line ( p : in Bracket_Polynomial );
  procedure put_line ( file : in file_type; p : in Bracket_Polynomial );

  -- DESCRIPTION :
  --   Writes the polynomial on standard output or on file,
  --   depending whether a file is an argument of the procedure,
  --   taking a new line for each term.

end QuadDobl_Bracket_Polynomials_io;
