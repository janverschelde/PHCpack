with text_io;                            use text_io;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Matrices;

generic

  with package Ring_io is new Abstract_Ring_io(<>);  use Ring_io.Ring;
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Matrices is new Generic_Matrices(Ring_io.Ring,Vectors);

package Generic_Polynomial_Systems_io is

-- DESCRIPTION :
--  This package contains routines for the input and output
--  of human-readable polynomial systems.
--  Note that for every write-procedure, there is a read-
--  procedure that can read the written system.
--  See also the package for the input and output
--  of polynomials in n unknowns.

-- THE INPUT OPERATIONS :

  procedure get ( n : in natural; s : out Poly_Sys );
  procedure get ( n,m : in natural; s : out Poly_Sys );
  procedure get ( file : in file_type; n : in natural; s : out Poly_Sys );
  procedure get ( file : in file_type; n,m : in natural; s : out Poly_Sys );
  procedure get ( s : out Poly_Sys );
  procedure get ( file : in file_type; s : out Poly_Sys );

  -- DESCRIPTION :
  --   A polynomial system is read; n polynomials are read.

  -- REQUIRED : 
  --  * all unknows must begin with a letter and may have
  --    no symbols like '+', '-', '*', '^', '/', ';' or brackets in them;
  --    i = sqrt(-1) is reserved for complex numbers representation
  --  * each symbol is limited to 3 characters
  --  * the input is terminated by the delimiter
  --  * no blanks may occur in the numbers
  --  * if a file is specified, then it must be opened for input
  --  * n and m should both be on the first line !
  --  * if n stands alone, it should be followed immediately by an 
  --    end-of-line symbol.
 
  -- NOTE :
  --   The end_of_line symbol is not read.

  -- ON ENTRY :
  --   file       file_type where the input is,
  --              if not specified, then standard input is assumed;
  --   n          the number of equations,
  --              if not specified, then n will first be read;
  --   m          the number of unknowns;
  --              if not specified, then m will first be read. 

  -- ON RETURN :
  --   s         a polynomial system.

-- MORE USER FRIENDLY INPUT OPERATIONS :

  procedure get ( lp : in out Link_to_Poly_Sys );
  procedure get ( file : in file_type; lp : in out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   The first routine asks for a file name, while the second one
  --   assumes everything is on file and nothing has to be read from
  --   standard input.

  -- NOTE :
  --   The end_of_line symbol is read at the end of the polynomial system.

-- THE OUTPUT OPERATIONS :

  procedure put ( n : out natural; s : in Poly_Sys; pow : in power := '*' );
  procedure put ( n,m : out natural; s : in Poly_Sys; pow : in power := '*' );
  procedure put ( file : in file_type;
                  n : out natural; s : in Poly_Sys; pow : in power := '*' );
  procedure put ( file : in file_type;
                  n,m : out natural; s : in Poly_Sys; pow : in power := '*' );
  procedure put ( s : in Poly_Sys; pow : in power );
  procedure put ( file : in file_type; s : in Poly_Sys; pow : in power );
  procedure put ( s : in Poly_Sys );
  procedure put ( file : in file_type; s : in Poly_Sys );

  -- DESCRIPTION :
  --   A polynomial system is written on standard output or on file.

  -- ON ENTRY :
  --   file       file where the output must come;
  --              if not specified, then standard output is assumed
  --   s          a polynomial system;
  --   pow        kind of power symbol used.

  -- ON RETURN :
  --   n          the number of equations of p,
  --              if not specified, n will first be written;
  --   m          the number of unknowns.
  --              if not specified, m will first be written.

  procedure put_line ( s : in Poly_Sys );
  procedure put_line ( file : in file_type; s : in Poly_Sys );
  procedure put_line ( s : in Poly_Sys; pow : in Power );
  procedure put_line ( file : in file_type; s : in Poly_Sys; pow : in Power );

  -- DESCRIPTION :
  --   Writes the polynomials, every term on a separate line.

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

end Generic_Polynomial_Systems_io;
