with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Standard_Complex_Poly_Systems_io is

-- DESCRIPTION :
--  This package contains routines for the i/o for polynomial systems.

-- THE INPUT OPERATIONS :

  procedure get ( n : in out natural32; p : in out Poly_Sys );
  procedure get ( n,m : in out natural32; p : in out Poly_Sys );
  procedure get ( file : in file_type;
                  n : in out natural32; p : in out Poly_Sys );
  procedure get ( file : in file_type;
                  n,m : in out natural32; p : in out Poly_Sys );
  procedure get ( p : in out Poly_Sys );
  procedure get ( file : in file_type; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   A polynomial system is read; n polynomials are read.

  -- ON ENTRY :
  --   file       file_type where the input is,
  --              if not specified, then standard input is assumed;
  --   n          the number of equations,
  --              if specified, then n will first be read;
  --   m          the number of unknowns;
  --              if specified, then m will first be read. 

  -- ON RETURN :
  --   p          a polynomial system.

-- MORE USER FRIENDLY INPUT OPERATIONS :

  procedure get ( lp : out Link_to_Poly_Sys );
  procedure get ( file : in file_type; lp : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   The first routine asks for a file name, while the second one
  --   assumes everything is on file and nothing has to be read from
  --   standard input.

  -- NOTE :
  --   The end_of_line symbol is read at the end of the polynomial system.

-- THE OUTPUT OPERATIONS :

  procedure put ( n : in natural32; p : in Poly_Sys; pow : in power := '^' );
  procedure put ( n,m : in natural32; p : in Poly_Sys; pow : in power := '^' );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly_Sys; pow : in power := '^' );
  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Poly_Sys; pow : in power := '^' );
  procedure put ( p : in Poly_Sys; pow : in power );
  procedure put ( file : in file_type; p : in Poly_Sys; pow : in power );
  procedure put ( p : in Poly_Sys );
  procedure put ( file : in file_type; p : in Poly_Sys );

  procedure put ( n : in natural32; p : in Poly_Sys;
                  s : in Array_of_Symbols; pow : in power := '^' );
  procedure put ( n,m : in natural32; p : in Poly_Sys;
                  s : in Array_of_Symbols; pow : in power := '^' );
  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly_Sys;
                  s : in Array_of_Symbols; pow : in power := '^' );
  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Poly_Sys;
                  s : in Array_of_Symbols; pow : in power := '^' );
  procedure put ( p : in Poly_Sys; s : in Array_of_Symbols; pow : in power );
  procedure put ( file : in file_type;
                  p : in Poly_Sys; s : in Array_of_Symbols; pow : in power );
  procedure put ( p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put ( file : in file_type;
                  p : in Poly_Sys; s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   A polynomial system is written on standard output or on file.

  -- ON ENTRY :
  --   file       file where the output must come,
  --              if not specified, then standard output is assumed;
  --   n          the number of equations of p,
  --              if specified, n will first be written;
  --   m          the number of unknowns,
  --              if specified, m will first be written;
  --   p          a polynomial system;
  --   s          symbols to be used for the unknowns;
  --   pow        kind of power symbol used.

  procedure put_line ( p : in Poly_Sys );
  procedure put_line ( file : in file_type; p : in Poly_Sys );
  procedure put_line ( p : in Poly_Sys; pow : in Power );
  procedure put_line ( file : in file_type; p : in Poly_Sys; pow : in Power );
  procedure put_line ( p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put_line ( file : in file_type;
                       p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put_line ( p : in Poly_Sys;
                       s : in Array_of_Symbols; pow : in Power );
  procedure put_line ( file : in file_type; p : in Poly_Sys; 
                       s : in Array_of_Symbols; pow : in Power );

  -- DESCRIPTION :
  --   Writes the polynomials, every term on a separate line,
  --   which is useful for random complex coefficient systems.

  procedure put_pair ( p : in Poly_Sys );
  procedure put_pair ( file : in file_type; p : in Poly_Sys );
  procedure put_pair ( p : in Poly_Sys; pow : in Power );
  procedure put_pair ( file : in file_type; p : in Poly_Sys; pow : in Power );
  procedure put_pair ( p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put_pair ( file : in file_type;
                       p : in Poly_Sys; s : in Array_of_Symbols );
  procedure put_pair ( p : in Poly_Sys;
                       s : in Array_of_Symbols; pow : in Power );
  procedure put_pair ( file : in file_type; p : in Poly_Sys; 
                       s : in Array_of_Symbols; pow : in Power );

  -- DESCRIPTION :
  --   Writes the polynomials with the terms in pairs on separate lines,
  --   which is useful when all coefficients are general double floats.

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

end Standard_Complex_Poly_Systems_io;
