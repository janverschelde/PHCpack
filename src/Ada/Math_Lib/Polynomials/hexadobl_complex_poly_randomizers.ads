with HexaDobl_Complex_Polynomials;       use HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Systems;      use HexaDobl_Complex_Poly_Systems;

package HexaDobl_Complex_Poly_Randomizers is

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

end HexaDobl_Complex_Poly_Randomizers;
