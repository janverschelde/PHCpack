with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with TripDobl_Complex_Laur_Systems;     use TripDobl_Complex_Laur_Systems;

package TripDobl_Complex_Laur_Systems_io is

-- DESCRIPTION :
--   This package provides basic input/output routines of Laurent polynomial
--   systems with triple double complex coefficients.

  procedure get ( p : out Link_to_Laur_Sys );
  procedure get ( file : in file_type; p : out Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   Reads in a system with multiprecision coefficients and converts
  --   the result to a system with triple double complex coefficients.

  procedure put ( p : in Laur_Sys );
  procedure put ( file : in file_type; p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file.   
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.

  procedure put ( n : in natural32; p : in Laur_Sys );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file,
  --   after first writing the number n of polynomials in p.
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.

  procedure put ( n,m : in natural32; p : in Laur_Sys );
  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file,
  --   after first writing the number n of polynomials in p,
  --   followed by the number m of variables in each polynomial.
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.

  procedure put_line ( p : in Laur_Sys );
  procedure put_line ( file : in file_type; p : in Laur_Sys );

  -- DESCRIPTION :
  --   A new line is started for each new monomial.

end TripDobl_Complex_Laur_Systems_io;
