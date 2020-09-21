with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;                      use Symbol_Table;
with OctoDobl_Complex_Poly_Systems;     use OctoDobl_Complex_Poly_Systems;

package OctoDobl_Complex_Poly_Systems_io is

-- DESCRIPTION :
--   This package provides basic input/output routines of polynomial systems
--   with octo double complex coefficients.

  procedure get ( p : out Link_to_Poly_Sys );
  procedure get ( file : in file_type; p : out Poly_Sys );
  procedure get ( file : in file_type; p : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Reads in a system with multiprecision coefficients and converts
  --   the result to a system with octo double complex coefficients.

  procedure put ( p : in Poly_Sys );
  procedure put ( file : in file_type; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file.   
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.

  procedure put ( n : in natural32; p : in Poly_Sys );
  procedure put ( file : in file_type; n : in natural32; p : in Poly_Sys );

  -- DESCRIPTION
  --   Writes n as the number of equation before writing p.

  procedure put ( n,m : in natural32; p : in Poly_Sys );
  procedure put ( file : in file_type; n,m : in natural32; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes the number n of polynomials followed by the number m 
  --   of variables on the same first line before writing p.

  procedure put_line ( p : in Poly_Sys );
  procedure put_line ( file : in file_type; p : in Poly_Sys );

  -- DESCRIPTION :
  --   A new line is started for each new monomial.

  procedure put ( p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  p : in Poly_Sys; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Writes the polynomial system p to standard output or to file,
  --   using the array of symbols as the symbols for the variables
  --   instead of the symbols stored in the symbol table.

end OctoDobl_Complex_Poly_Systems_io;
