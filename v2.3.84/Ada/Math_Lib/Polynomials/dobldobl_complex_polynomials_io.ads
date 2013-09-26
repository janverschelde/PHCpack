with text_io;                            use text_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;

package DoblDobl_Complex_Polynomials_io is

-- DESCRIPTION :
--   This package provides very basic output routines for polynomials
--   with double double complex coefficients.

  procedure put ( p : in Poly );
  procedure put ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   A polynomial in n unknowns is written on file or on standard output.

  -- ON ENTRY :
  --   file       file where the output must come,
  --              if not specified, then standard output is assumed;
  --   p          a polynomial in n unknows.

  procedure put_line ( p : in Poly );
  procedure put_line ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   Every term of the polynomial p will be written on a separate line.
  --   This is useful for polynomials with random complex coefficients.

end DoblDobl_Complex_Polynomials_io;
