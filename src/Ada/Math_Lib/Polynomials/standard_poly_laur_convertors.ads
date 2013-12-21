with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Standard_Poly_Laur_Convertors is

-- DESCRIPTION :
--   This package contains routines for converting ordinary polynomials
--   into Laurent polynomials.

  function Polynomial_to_Laurent_Polynomial
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Transforms a polynomial into a Laurent polynomial.

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Transforms a polynomial system into a Laurent polynomial system.

end Standard_Poly_Laur_Convertors;
