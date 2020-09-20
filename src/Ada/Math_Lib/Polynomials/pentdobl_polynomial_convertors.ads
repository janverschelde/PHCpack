with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Penta_Double_Polynomials;
with Penta_Double_Poly_Systems;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Laurentials;
with PentDobl_Complex_Poly_Systems;
with PentDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package PentDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with penta double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Penta_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Penta_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Penta_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Penta_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Penta_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Penta_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Penta_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Penta_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with penta double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_PentDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return PentDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_PentDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_PentDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return PentDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_PentDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex penta double coefficients.

  function Standard_Laurential_to_PentDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return PentDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_PentDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return PentDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_PentDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return PentDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_PentDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return PentDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex penta double coefficients.

  function Penta_Double_to_Standard_Polynomial
             ( p : Penta_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Penta_Double_to_Standard_Poly_Sys
             ( p : Penta_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Penta_Double_to_Multprec_Polynomial
             ( p : Penta_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Penta_Double_to_Multprec_Poly_Sys
             ( p : Penta_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with penta double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the penta double coefficients matter.

  function PentDobl_Complex_to_Standard_Polynomial
             ( p : PentDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function PentDobl_Complex_to_Standard_Poly_Sys
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function PentDobl_Complex_to_Multprec_Polynomial
             ( p : PentDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function PentDobl_Complex_to_Multprec_Poly_Sys
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex penta double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the penta double coefficients matter.

  function PentDobl_Complex_to_Standard_Laurential
             ( p : PentDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function PentDobl_Complex_to_Standard_Laur_Sys
             ( p : PentDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function PentDobl_Complex_to_Multprec_Laurential
             ( p : PentDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function PentDobl_Complex_to_Multprec_Laur_Sys
             ( p : PentDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex penta double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the penta double coefficients matter.

end PentDobl_Polynomial_Convertors;
