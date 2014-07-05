with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;      use DoblDobl_Complex_Laur_Systems;

package DoblDobl_Poly_Laur_Convertors is

-- DESCRIPTION :
--   This package contains routines for converting ordinary polynomials
--   into Laurent polynomials.

  function Polynomial_to_Laurent_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Transforms a polynomial into a Laurent polynomial.

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Transforms a polynomial system into a Laurent polynomial system.

end DoblDobl_Poly_Laur_Convertors;
