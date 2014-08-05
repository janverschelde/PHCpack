with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Jaco_Matrices;    use DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Matrices;    use DoblDobl_Complex_Poly_Matrices;

package DoblDobl_Embed_Polynomials is

-- DESCRIPTION :
--   The operations in this package extend the polynomials
--   with extra variables.  This corresponds to the embedding
--   of the polynomials into a larger space.

  function Add_Variables ( p : Poly; k : natural32 ) return Poly;
  function Add_Variables ( p : Poly_Sys; k : natural32 ) return Poly_Sys;
  function Add_Variables ( p : Jaco_Mat; k : natural32 ) return Jaco_Mat;
  function Add_Variables ( p : Matrix; k : natural32 ) return Matrix;

  -- DESCRIPTION :
  --   Extends the representation of the degrees with k more variables,
  --   whose exponent equals 0.

end DoblDobl_Embed_Polynomials;
