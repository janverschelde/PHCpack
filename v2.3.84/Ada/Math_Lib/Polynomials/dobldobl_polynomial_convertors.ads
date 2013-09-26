with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Double_Double_Polynomials;
with Double_Double_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package DoblDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with double double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Double_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Double_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Double_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Double_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with double double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_DoblDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_DoblDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_DoblDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_DoblDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex double double coefficients.

  function Standard_Laurential_to_DoblDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return DoblDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_DoblDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_DoblDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return DoblDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_DoblDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex double double coefficients.

  function Double_Double_to_Standard_Polynomial
             ( p : Double_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Double_Double_to_Standard_Poly_Sys
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Double_Double_to_Multprec_Polynomial
             ( p : Double_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Double_Double_to_Multprec_Poly_Sys
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with double double coefficients into
  --   a polynomial with standard/multprec complex coefficients.
  --   Low parts of the double double coefficients are discarded
  --   when converting to standard coefficients.

  function DoblDobl_Complex_to_Standard_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function DoblDobl_Complex_to_Standard_Poly_Sys
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function DoblDobl_Complex_to_Multprec_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function DoblDobl_Complex_to_Multprec_Poly_Sys
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex double double coefficients
  --   into a polynomial with standard/multprec complex coefficients.
  --   Low parts of the double double coefficients are discarded
  --   when converting to standard coefficients.

  function DoblDobl_Complex_to_Standard_Laurential
             ( p : DoblDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function DoblDobl_Complex_to_Standard_Laur_Sys
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function DoblDobl_Complex_to_Multprec_Laurential
             ( p : DoblDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function DoblDobl_Complex_to_Multprec_Laur_Sys
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex double double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   Low parts of the double double coefficients are discarded
  --   when converting to standard coefficients.

end DoblDobl_Polynomial_Convertors;
