with text_io;                            use text_io;
with TripDobl_Complex_Laurentials;       use TripDobl_Complex_Laurentials;

package TripDobl_Complex_Laurentials_io is

-- DESCRIPTION :
--   This package provides very basic output routines for Laurent polynomials
--   with triple double complex coefficients.

-- THE INPUT OPERATIONS :

  procedure get ( p : out Poly );
  procedure get ( file : in file_type; p : out Poly );

  -- DESCRIPTION :
  --   Reads a multivariate polynomial from standard input or from file.

  -- ON ENTRY :
  --   file     optional file, must be opened for input.

  -- ON RETURN :
  --   p        Laurent polynomial in several variables.

-- THE OUTPUT OPERATIONS :

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

end TripDobl_Complex_Laurentials_io;
