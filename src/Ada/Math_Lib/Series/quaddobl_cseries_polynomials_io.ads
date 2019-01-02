with text_io;                           use text_io;
with QuadDobl_CSeries_Polynomials;      use QuadDobl_CSeries_Polynomials;

package QuadDobl_CSeries_Polynomials_io is

-- DESCRIPTION :
--   Provides basic output procedures for polynomials in several variables,
--   with series coefficients in quad double precision.

  procedure put ( t : in Term );
  procedure put ( file : in file_type; t : in Term );

  -- DESCRIPTION :
  --   Writes the term t to standard output or to file.
  --   The symbol 't' denotes the variable in the coefficient series,
  --   variables in the monomials are written with x1, x2, etc.

  procedure put ( p : in Poly );
  procedure put ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   Write the polynomial p to standard output or to file.

end QuadDobl_CSeries_Polynomials_io;
