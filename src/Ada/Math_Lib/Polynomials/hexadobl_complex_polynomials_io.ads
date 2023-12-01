with text_io;                            use text_io;
with Symbol_Table;                       use Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with HexaDobl_Complex_Numbers;           use HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Polynomials;       use HexaDobl_Complex_Polynomials;

package HexaDobl_Complex_Polynomials_io is

-- DESCRIPTION :
--   This package provides very basic output routines for polynomials
--   with hexa double complex coefficients.

  procedure get ( p : out Poly );
  procedure get ( file : in file_type; p : in out Poly );

  -- DESCRIPTION :
  --   Reads a multivariate polynomial from standard input or from file.

  procedure Write_Number ( file : in file_type; c : in Complex_Number );

  -- DESCRIPTION :
  --   Writes the hexa double complex number to file
  --   in the format ( real_part(c) + imag_part(c)*I ).

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

  procedure put_line ( p : in Poly; s : in Array_of_Symbols );
  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Every term of the polynomial p will be written on a separate line.
  --   This is useful for polynomials with random complex coefficients.
  --   Instead of the symbol table, the array of symbols will be used
  --   for the names of the variables.

  -- REQUIRED : s'range = 1..Number_of_Unknowns(p).

  procedure put ( p : in Poly; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  p : in Poly; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Writes the polynomial p to standard output or to file,
  --   using the array of symbols as the symbols for the variables
  --   instead of the symbols stored in the symbol table.

  -- REQUIRED : s'range = 1..Number_of_Unknowns(p).

end HexaDobl_Complex_Polynomials_io;
