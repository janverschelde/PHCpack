with text_io;                           use text_io;
with QuadDobl_Complex_Monomials;        use QuadDobl_Complex_Monomials;

package QuadDobl_Complex_Monomials_io is

-- DESCRIPTION :
--   A basic output procedure is defined by this package.

  procedure put ( m : in Monomial );
  procedure put ( file : in file_type; m : in Monomial );

  -- DESCRIPTION :
  --   Writes the contents of the monomial.

end QuadDobl_Complex_Monomials_io;
