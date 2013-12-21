with Multprec_Floating_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Floating_Poly_Systems;

package Multprec_Complex_to_Real_Poly is

-- DESCRIPTION :
--   This package offers type conversions between polynomials with
--   multiprecision complex coefficients into real ones, as needed
--   because the input/output is primarily for complex coefficients.

  function Convert_Complex_to_Real
             ( p : Multprec_Complex_Polynomials.Poly )
             return Multprec_Floating_Polynomials.Poly;
  function Convert_Complex_to_Real
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Floating_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Retains only the real part of every coefficient,
  --   assuming all coefficients are real, i.e.: after dropping the
  --   imaginary part of the complex coefficient nothing is lost.

  function Convert_Real_to_Complex
             ( p : Multprec_Floating_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Convert_Real_to_Complex
             ( p : Multprec_Floating_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   To every real coefficient a zero imaginary part is added
  --   so the polynomial (system) on return is of the right type.

end Multprec_Complex_to_Real_Poly;
