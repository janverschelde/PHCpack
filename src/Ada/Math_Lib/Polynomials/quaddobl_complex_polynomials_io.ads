with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;

package QuadDobl_Complex_Polynomials_io is

-- DESCRIPTION :
--   This package provides very basic output routines for polynomials
--   with quad double complex coefficients.

  procedure get ( p : in out Poly );
  procedure get ( file : in file_type; p : in out Poly );

  -- DESCRIPTION :
  --   Reads a multivariate polynomial from standard input or from file.

  procedure put ( p : in Poly );
  procedure put ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   A polynomial in n unknowns is written on file or on standard output.

  -- ON ENTRY :
  --   file       file where the output must come,
  --              if not specified, then standard output is assumed;
  --   p          a polynomial in n unknows.

  procedure put ( p : in Poly; dp : in natural32 ); 
  procedure put ( file : in file_type; p : in Poly; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the coefficients of the polynomial with as many decimal places
  --   as the value of dp, as needed to instantiate the abstract_ring_io.

  procedure put_line ( p : in Poly );
  procedure put_line ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   Every term of the polynomial p will be written on a separate line.
  --   This is useful for polynomials with random complex coefficients.

end QuadDobl_Complex_Polynomials_io;
