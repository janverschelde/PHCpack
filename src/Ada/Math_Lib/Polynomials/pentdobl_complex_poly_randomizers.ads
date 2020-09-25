with PentDobl_Complex_Polynomials;       use PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Systems;      use PentDobl_Complex_Poly_Systems;

package PentDobl_Complex_Poly_Randomizers is

-- DESCRIPTION :
--   This package offers routines for randomizing and perturbing
--   the coefficients of polynomials.
--   Except for the last three functions, the monomial structure 
--   remains the same, only random (real or complex) coefficients 
--   will replace the existing ones.

  function Complex_Randomize ( p : Poly ) return Poly;
  function Complex_Randomize ( p : Poly_Sys ) return Poly_Sys;
 
  -- DESCRIPTION :
  --   The real and imaginary part of the randomly generated
  --   coefficients are in [-1.0,1.0]

  function Complex_Randomize1 ( p : Poly ) return Poly;
  function Complex_Randomize1 ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Generates random complex coefficients with modulus one.

end PentDobl_Complex_Poly_Randomizers;
