with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Deca_Double_Polynomials;
with Deca_Double_Poly_Systems;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Laurentials;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package DecaDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with deca double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Deca_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Deca_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Deca_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Deca_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Deca_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Deca_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Deca_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Deca_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with deca double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_DecaDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return DecaDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_DecaDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_DecaDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return DecaDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_DecaDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex deca double coefficients.

  function Standard_Laurential_to_DecaDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return DecaDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_DecaDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return DecaDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_DecaDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return DecaDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_DecaDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return DecaDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex deca double coefficients.

  function Deca_Double_to_Standard_Polynomial
             ( p : Deca_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Deca_Double_to_Standard_Poly_Sys
             ( p : Deca_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Deca_Double_to_Multprec_Polynomial
             ( p : Deca_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Deca_Double_to_Multprec_Poly_Sys
             ( p : Deca_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with deca double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the deca double coefficients matter.

  function DecaDobl_Complex_to_Standard_Polynomial
             ( p : DecaDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function DecaDobl_Complex_to_Standard_Poly_Sys
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function DecaDobl_Complex_to_Multprec_Polynomial
             ( p : DecaDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function DecaDobl_Complex_to_Multprec_Poly_Sys
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex deca double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the deca double coefficients matter.

  function DecaDobl_Complex_to_Standard_Laurential
             ( p : DecaDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function DecaDobl_Complex_to_Standard_Laur_Sys
             ( p : DecaDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function DecaDobl_Complex_to_Multprec_Laurential
             ( p : DecaDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function DecaDobl_Complex_to_Multprec_Laur_Sys
             ( p : DecaDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex deca double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the deca double coefficients matter.

end DecaDobl_Polynomial_Convertors;
