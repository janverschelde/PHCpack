with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Hexa_Double_Polynomials;
with Hexa_Double_Poly_Systems;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Laurentials;
with HexaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package HexaDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with hexa double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Hexa_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Hexa_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Hexa_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Hexa_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Hexa_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Hexa_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Hexa_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Hexa_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with hexa double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_HexaDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_HexaDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_HexaDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_HexaDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex hexa double coefficients.

  function Standard_Laurential_to_HexaDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return HexaDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_HexaDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return HexaDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_HexaDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return HexaDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_HexaDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return HexaDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex hexa double coefficients.

  function Hexa_Double_to_Standard_Polynomial
             ( p : Hexa_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Hexa_Double_to_Standard_Poly_Sys
             ( p : Hexa_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Hexa_Double_to_Multprec_Polynomial
             ( p : Hexa_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Hexa_Double_to_Multprec_Poly_Sys
             ( p : Hexa_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with hexa double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the hexa double coefficients matter.

  function HexaDobl_Complex_to_Standard_Polynomial
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function HexaDobl_Complex_to_Standard_Poly_Sys
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function HexaDobl_Complex_to_Multprec_Polynomial
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function HexaDobl_Complex_to_Multprec_Poly_Sys
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex hexa double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the hexa double coefficients matter.

  function HexaDobl_Complex_to_Standard_Laurential
             ( p : HexaDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function HexaDobl_Complex_to_Standard_Laur_Sys
             ( p : HexaDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function HexaDobl_Complex_to_Multprec_Laurential
             ( p : HexaDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function HexaDobl_Complex_to_Multprec_Laur_Sys
             ( p : HexaDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex hexa double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the hexa double coefficients matter.

end HexaDobl_Polynomial_Convertors;
