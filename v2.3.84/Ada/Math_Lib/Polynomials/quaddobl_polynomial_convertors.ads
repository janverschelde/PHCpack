with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Quad_Double_Polynomials;
with Quad_Double_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package QuadDobl_Polynomial_Convertors is

-- DESCRIPTION :
--   Polynomials with quad double coefficients can be converted from
--   and into polynomials with standard (hardware) coefficients,
--   or into polynomials with arbitrary multiprecision coefficients.

  function Standard_Polynomial_to_Double_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Quad_Double_Polynomials.Poly;
  function Standard_Poly_Sys_to_Double_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Quad_Double_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_Double_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Quad_Double_Polynomials.Poly;
  function Multprec_Poly_Sys_to_Double_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Quad_Double_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Turns a polynomial with standard/multprec complex coefficients
  --   into a polynomial with quad double coefficients.
  --   Imaginary parts of the coefficients of p are discarded.

  function Standard_Polynomial_to_QuadDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Standard_Poly_Sys_to_QuadDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Polynomial_to_QuadDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Multprec_Poly_Sys_to_QuadDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Turns a polynomial with standard/multprec complex coefficients
  --   into a polynomial with complex quad double coefficients.

  function Standard_Laurential_to_QuadDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly;
  function Standard_Laur_Sys_to_QuadDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function Multprec_Laurential_to_QuadDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly;
  function Multprec_Laur_Sys_to_QuadDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Turns a Laurent polynomial with standard/multprec complex coefficients
  --   into a Laurent polynomial with complex quad double coefficients.

  function Quad_Double_to_Standard_Polynomial
             ( p : Quad_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function Quad_Double_to_Standard_Poly_Sys
             ( p : Quad_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Quad_Double_to_Multprec_Polynomial
             ( p : Quad_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function Quad_Double_to_Multprec_Poly_Sys
             ( p : Quad_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with quad double coefficients into
  --   a polynomial with standard complex coefficients.
  --   Low parts of the quad double coefficients are discarded
  --   when converting to standard coefficients.

  function QuadDobl_Complex_to_Standard_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;
  function QuadDobl_Complex_to_Standard_Poly_Sys
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function QuadDobl_Complex_to_Multprec_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly;
  function QuadDobl_Complex_to_Multprec_Poly_Sys
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Converts a polynomial with complex quad double coefficients
  --   into a polynomial with standard complex coefficients.
  --   Low parts of the quad double coefficients are discarded
  --   when converting to standard coefficients.

  function QuadDobl_Complex_to_Standard_Laurential
             ( p : QuadDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly;
  function QuadDobl_Complex_to_Standard_Laur_Sys
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function QuadDobl_Complex_to_Multprec_Laurential
             ( p : QuadDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly;
  function QuadDobl_Complex_to_Multprec_Laur_Sys
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial with complex quad double coefficients
  --   into a Laurent polynomial with standard complex coefficients.
  --   Low parts of the quad double coefficients are discarded
  --   when converting to standard coefficients.

end QuadDobl_Polynomial_Convertors;
