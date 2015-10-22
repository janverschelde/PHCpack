with Standard_Bracket_Polynomials;
with DoblDobl_Bracket_Polynomials;
with QuadDobl_Bracket_Polynomials;

package Bracket_Polynomial_Convertors is

-- DESCRIPTION :
--   Converts bracket polynomials with coefficients in standard double
--   precision into brackt polynomials with coefficients in double double
--   or quad double precision.

  function Convert
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial )
             return DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
  function Convert
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial )
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial;

  -- DESCRIPTION :
  --   On return is a bracket polynomial with the same monomials,
  --   but with coefficients in double double or quad double precision.

end Bracket_Polynomial_Convertors;
