with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;

package Standard_Complex_Prod_Systems_io is

-- DESCRIPTION :
--   Input and output routines for systems of products of polynomials.

  procedure get ( lp : out Link_to_Prod_Sys );
  procedure get ( file : in file_type; lp : out Link_to_Prod_Sys );

  -- DESCRIPTION :
  --   Asks the user if the system is available on file, and read then
  --   the product system either from file, or from standard input.
  --   The symbol table will be initialized if necessary.

  procedure get ( p : out Prod_Sys );
  procedure get ( file : in file_type; p : out Prod_Sys );

  -- DESCRIPTION :
  --   Reads a system of products of polynomials from standard input,
  --   or from file, following the format of product polynomials.
  --   Every polynomial in the product must be enclosed between
  --   round brackets with a '*' in between two factors.
  --   The symbol terminating the product is a semicolon.
 
  -- REQUIRED :
  --   The symbol table must be initialized accordingly.

  procedure put ( p : in Prod_Sys );
  procedure put ( n : in natural32; p : in Prod_Sys );
  procedure put ( n,m : in natural32; p : in Prod_Sys );
  procedure put ( file : in file_type; p : in Prod_Sys );
  procedure put ( file : in file_type; n : in natural32; p : in Prod_Sys );
  procedure put ( file : in file_type; n,m : in natural32; p : in Prod_Sys );

  -- DESCRIPTION :
  --   Writes a product polynomial to screen or file,
  --   in the same format as described above with the get.
  --   If provided as n and m, then the number of equations n
  --   and the number of unknowns will be written to file.

  procedure put_line ( p : in Prod_Sys );
  procedure put_line ( n : in natural32; p : in Prod_Sys );
  procedure put_line ( n,m : in natural32; p : in Prod_Sys );
  procedure put_line ( file : in file_type; p : in Prod_Sys );
  procedure put_line ( file : in file_type;
                       n : in natural32; p : in Prod_Sys );
  procedure put_line ( file : in file_type;
                       n,m : in natural32; p : in Prod_Sys );

  -- DESCRIPTION :
  --   Writes a product polynomial to screen or file,
  --   using a new line for every term in the product polynomials.

end Standard_Complex_Prod_Systems_io;
