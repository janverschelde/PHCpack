with text_io;                           use text_io;
with DoblDobl_Complex_Laur_Systems;     use DoblDobl_Complex_Laur_Systems;

package DoblDobl_Complex_Laur_Systems_io is

-- DESCRIPTION :
--   This package provides basic input/output routines of Laurent polynomial
--   systems with double double complex coefficients.

  procedure get ( p : out Link_to_Laur_Sys );
  procedure get ( file : in file_type; p : out Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   Reads in a system with multiprecision coefficients and converts
  --   the result to a system with double double complex coefficients.

  procedure put ( p : in Laur_Sys );
  procedure put ( file : in file_type; p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file.   
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.

  procedure put_line ( p : in Laur_Sys );
  procedure put_line ( file : in file_type; p : in Laur_Sys );

  -- DESCRIPTION :
  --   A new line is started for each new monomial.

end DoblDobl_Complex_Laur_Systems_io;
