with Double_Double_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with Double_Double_Poly_Systems;

package DoblDobl_Complex_to_Real_Poly is

-- DESCRIPTION :
--   This package offers type conversions between polynomials with
--   complex coefficients into real ones, as needed because
--   the input/output is defined primarily for complex coefficients.

  function Is_Real ( p : DoblDobl_Complex_Polynomials.Poly ) return boolean;
  function Is_Real
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all coefficients in p have a zero imaginary part.

  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly;
  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys;
  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
             return Double_Double_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Retains only the real part of every coefficient,
  --   assuming all coefficients are real, i.e.: after dropping the
  --   imaginary part of the complex coefficient nothing is lost.

  function Convert_Real_to_Complex
             ( p : Double_Double_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly;
  function Convert_Real_to_Complex
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Convert_Real_to_Complex
             ( p : Double_Double_Poly_Systems.Link_to_Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   To every real coefficient a zero imaginary part is added
  --   so the polynomial (system) on return is of the right type.

end DoblDobl_Complex_to_Real_Poly;
