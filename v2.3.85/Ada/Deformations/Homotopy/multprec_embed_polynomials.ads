with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;     use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Jaco_Matrices;    use Multprec_Complex_Jaco_Matrices;

package Multprec_Embed_Polynomials is

-- DESCRIPTION :
--   The operations in this package extend the polynomials
--   with extra variables.  This corresponds to the embedding
--   of the polynomials into a larger space.

  function Add_Variables ( p : Poly; k : natural32 ) return Poly;
  function Add_Variables ( p : Poly_Sys; k : natural32 ) return Poly_Sys;
  function Add_Variables ( p : Jaco_Mat; k : natural32 ) return Jaco_Mat;

  -- DESCRIPTION :
  --   Extends the representation of the degrees with k more variables,
  --   whose exponent equals 0.

end Multprec_Embed_Polynomials;
