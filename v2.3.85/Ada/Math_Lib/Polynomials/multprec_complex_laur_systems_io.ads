with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;                      use Symbol_Table;
with Multprec_Complex_Laur_Systems;     use Multprec_Complex_Laur_Systems;

package Multprec_Complex_Laur_Systems_io is

-- DESCRIPTION :
--   This package provides input/output routines of Laurent polynomial systems
--   with multiprecision complex coefficients.

-- INPUT ROUTINES :

  procedure get ( p : in out Laur_Sys );
  procedure get ( file : in file_type; p : in out Laur_Sys );

  -- DESCRIPTION :
  --   Will read as many polynomials from standard input or from file
  --   as the range of p.

  procedure get ( n : in out natural32; p : in out Laur_Sys ); 
  procedure get ( file : in file_type;
                  n : in out natural32; p : in out Laur_Sys ); 

  -- DESCRIPTION : 
  --   First reads the number of polynomials: n from standard input or from
  --   file.  Then the reading of the polynomials starts, up to n or until
  --   the range of p is full.

  procedure get ( lp : out Link_to_Laur_Sys );
  procedure get ( file : in file_type; lp : out Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   The first routine asks for a file name, while the second one
  --   assumes everything is on file and nothing has to be read from
  --   standard input.

  -- NOTE :
  --   The end_of_line symbol is read at the end of the polynomial system.

-- OUTPUT ROUTINES

  procedure put ( p : in Laur_Sys );
  procedure put ( n : in natural32; p : in Laur_Sys );
  procedure put ( file : in file_type; p : in Laur_Sys );
  procedure put ( file : in file_type; n : in natural32; p : in Laur_Sys );
  procedure put ( p : in Laur_Sys; s : in Array_of_Symbols );
  procedure put ( n : in natural32; p : in Laur_Sys; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  p : in Laur_Sys; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Laur_Sys; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Writes the polynomials to screen or file.   
  --   The polynomials are separated by semicolons and a new line
  --   is started for each new polynomial.
  --   If n is provided, then the value for n will be first written.

  procedure put_line ( p : in Laur_Sys );
  procedure put_line ( file : in file_type; p : in Laur_Sys );
  procedure put_line ( p : in Laur_Sys; s : in Array_of_Symbols );
  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   A new line is started for each new monomial.

end Multprec_Complex_Laur_Systems_io;
