with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;

package Multprec_Complex_Polynomials_io is

-- DESCRIPTION :
--   This package contains routines for the input and output of
--   multi-precision complex multivariate polynomials in symbolic form.

-- EXAMPLES : i denotes sqrt(-1)
--   ex1 : x**2*y + 1/2*z*y**2 - 2*z + y**3 + x - 1E9/-8.E-6* y + 3;
--   ex2 : x^2*y + z*y^2 - 2*z + y^3 + x - y + 3;
--   ex3 : (1.01 + 2.8*i)*x1**2*x2 + x3**2*x1 - 3*x1 + 2*x2*x3 - 3;
--   ex4 : (x1^2*x2 + x3^2*x1 - 3*x1 + 2*x2*x3 - 3)*x2**2*(x2-1+i);

-- SETTING THE WORKING PRECISION :

  procedure Set_Working_Precision ( size : in natural32 );

  -- DESCRIPTION :
  --   Sets the working precision to evaluate fractions in the get
  --   procedures below.

-- THE INPUT OPERATIONS :

  procedure get ( n : in out natural32; p : in out Poly );
  procedure get ( file : in file_type; n : in out natural32; p : in out Poly );
  procedure get ( p : in out Poly );
  procedure get ( file : in file_type; p : in out Poly );

  -- DESCRIPTION :
  --   A polynomial will be read from file or standard input.
  --   No ordening on the monomials is assumed.

  -- REQUIRED : 
  --  * all unknows must begin with a letter and may have
  --    no symbols like '+', '-', '*', '^', '/', ';' or brackets in them,
  --    i = sqrt(-1) is reserved for complex numbers representation;
  --  * each symbol is limited to 3 characters;
  --  * the input is terminated by the delimiter;
  --  * no blanks may occur in the numbers;
  --  * if specified, then the file must already be opened for input.

  -- NOTE :
  --   The end_of_line symbol is not read.

  -- ON ENTRY :
  --   file       file where the input is,
  --              if not specified, then standard input is assumed;
  --   n          the number of unknowns,
  --              if specified, then n will first be read, otherwise it is
  --              is assumed that the symbol table is already initialized.

  -- ON RETURN :
  --   p          a polynomial in n unknowns.

-- THE OUTPUT OPERATIONS :

  procedure put ( p : in Poly; pow : in Power );
  procedure put ( file : in file_type; p : in Poly; pow : in Power );
  procedure put ( n : in natural32; p : in Poly; pow : in Power );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly; pow : in Power );
  procedure put ( p : in Poly );
  procedure put ( file : in file_type; p : in Poly );

  procedure put ( p : in Poly; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  p : in Poly; s : in Array_of_Symbols );
  procedure put ( n : in natural32; p : in Poly; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   A polynomial in n unknowns is written on file or on standard output.

  -- REQUIRED :
  --   If specified, then the file is must already been opened for output.
  --   If provided, s'length >= Number_of_Unknowns(p).

  -- ON ENTRY :
  --   file       file where the output must come,
  --              if not specified, then standard output is assumed;
  --   p          a polynomial in n unknowns;
  --   s          symbols to be used for the unknowns;
  --   pow        kind of power symbol used, the default is '*'.

  -- ON RETURN :
  --   n          the number of unknows of p,
  --              if specified, n will first be written.

  procedure put_terms ( p : in Poly; dp : in natural32 );
  procedure put_terms ( file : in file_type; p : in Poly; dp : in natural32 );
  procedure put ( p : in Poly; dp : in natural32 );
  procedure put ( file : in file_type; p : in Poly; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes the coefficients of the polynomials with dp decimal places.
  --   The "put_terms" procedure does not write the delimiter at the end.

  procedure put_line ( p : in Poly );
  procedure put_line ( file : in file_type; p : in Poly );
  procedure put_line ( p : in Poly; pow : in Power );
  procedure put_line ( file : in file_type; p : in Poly; pow : in Power );
  procedure put_line ( p : in Poly; s : in Array_of_Symbols );
  procedure put_line ( file : in file_type;
                       p : in Poly; s : in Array_of_Symbols );
  procedure put_line ( p : in Poly; s : in Array_of_Symbols; pow : in Power );
  procedure put_line ( file : in file_type;
                       p : in Poly; s : in Array_of_Symbols; pow : in Power );

  -- DESCRIPTION :
  --   Every term of the polynomial p will be written on a separate line.
  --   This is useful for polynomials with random complex coefficients.

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

end Multprec_Complex_Polynomials_io;
