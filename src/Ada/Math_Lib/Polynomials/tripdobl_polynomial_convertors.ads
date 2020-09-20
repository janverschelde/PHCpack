with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Triple_Double_Polynomials;
with Triple_Double_Poly_Systems;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Laurentials;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package TripDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with triple double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Triple_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Triple_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Triple_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Triple_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Triple_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Triple_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Triple_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Triple_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with triple double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_TripDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return TripDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_TripDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_TripDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return TripDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_TripDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex triple double coefficients.

  function Standard_Laurential_to_TripDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return TripDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_TripDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return TripDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_TripDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return TripDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_TripDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return TripDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex triple double coefficients.

  function Triple_Double_to_Standard_Polynomial
             ( p : Triple_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Triple_Double_to_Standard_Poly_Sys
             ( p : Triple_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Triple_Double_to_Multprec_Polynomial
             ( p : Triple_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Triple_Double_to_Multprec_Poly_Sys
             ( p : Triple_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with triple double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the triple double coefficients matter.

  function TripDobl_Complex_to_Standard_Polynomial
             ( p : TripDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function TripDobl_Complex_to_Standard_Poly_Sys
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function TripDobl_Complex_to_Multprec_Polynomial
             ( p : TripDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function TripDobl_Complex_to_Multprec_Poly_Sys
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex triple double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the triple double coefficients matter.

  function TripDobl_Complex_to_Standard_Laurential
             ( p : TripDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function TripDobl_Complex_to_Standard_Laur_Sys
             ( p : TripDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function TripDobl_Complex_to_Multprec_Laurential
             ( p : TripDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function TripDobl_Complex_to_Multprec_Laur_Sys
             ( p : TripDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex triple double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   When converting to standard coefficients,
  --   only the highest part of the triple double coefficients matter.

end TripDobl_Polynomial_Convertors;
