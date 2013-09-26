with text_io;                            use text_io;
with Standard_Complex_Poly_Lists;        use Standard_Complex_Poly_Lists;

package Standard_Complex_Poly_Lists_io is

-- DESCRIPTION :
--   Input and output routines of products of polynomials.

  procedure get ( p : out Prod_Poly );
  procedure get ( file : in file_type; p : out Prod_Poly );

  -- DESCRIPTION :
  --   Reads a product of polynomials from standard input or from file.
  --   Every polynomial in the product must be enclosed between
  --   round brackets with a '*' in between two factors.
  --   The symbol terminating the product is a semicolon.

  -- REQUIRED :
  --   The symbol table is initialized with the number of symbols
  --   expected to appear in the list of polynomials.

  procedure put ( p : in Prod_Poly );
  procedure put ( file : in file_type; p : in Prod_Poly );

  -- DESCRIPTION :
  --   Writes the product polynomial to screen or to file,
  --   in the same format as described above with the get.

  procedure put_line ( p : in Prod_Poly );
  procedure put_line ( file : in file_type; p : in Prod_Poly );

  -- DESCRIPTION :
  --   Writes the product polynomial to screen or to file,
  --   using one new line for every term in every polynomial.

end Standard_Complex_Poly_Lists_io;
