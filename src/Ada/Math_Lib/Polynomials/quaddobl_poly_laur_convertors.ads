with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;      use QuadDobl_Complex_Laur_Systems;

package QuadDobl_Poly_Laur_Convertors is

-- DESCRIPTION :
--   This package contains routines for converting ordinary polynomials
--   into Laurent polynomials.

  function Polynomial_to_Laurent_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Transforms a polynomial into a Laurent polynomial.

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Transforms a polynomial system into a Laurent polynomial system.

end QuadDobl_Poly_Laur_Convertors;
