with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Octo_Double_Polynomials;
with Octo_Double_Poly_Systems;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Laurentials;
with OctoDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package OctoDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with octo double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Octo_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Octo_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Octo_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Octo_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Octo_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Octo_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Octo_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Octo_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with octo double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_OctoDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return OctoDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_OctoDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_OctoDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return OctoDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_OctoDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex octo double coefficients.

  function Standard_Laurential_to_OctoDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return OctoDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_OctoDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return OctoDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_OctoDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return OctoDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_OctoDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return OctoDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex octo double coefficients.

  function Octo_Double_to_Standard_Polynomial
             ( p : Octo_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Octo_Double_to_Standard_Poly_Sys
             ( p : Octo_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Octo_Double_to_Multprec_Polynomial
             ( p : Octo_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Octo_Double_to_Multprec_Poly_Sys
             ( p : Octo_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with octo double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the octo double coefficients matter.

  function OctoDobl_Complex_to_Standard_Polynomial
             ( p : OctoDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function OctoDobl_Complex_to_Standard_Poly_Sys
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function OctoDobl_Complex_to_Multprec_Polynomial
             ( p : OctoDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function OctoDobl_Complex_to_Multprec_Poly_Sys
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex octo double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the octo double coefficients matter.

  function OctoDobl_Complex_to_Standard_Laurential
             ( p : OctoDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function OctoDobl_Complex_to_Standard_Laur_Sys
             ( p : OctoDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function OctoDobl_Complex_to_Multprec_Laurential
             ( p : OctoDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function OctoDobl_Complex_to_Multprec_Laur_Sys
             ( p : OctoDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex octo double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the octo double coefficients matter.

end OctoDobl_Polynomial_Convertors;
