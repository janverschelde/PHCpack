with Standard_Complex_Polynomials;
with Standard_Series_Polynomials;

package Series_and_Polynomials is

-- DESCRIPTION :
--   Exports functions to convert polynomials with complex coefficients
--   into polynomials with series as coefficients, and vice versa.
--   The conversion routines give immediate access to symbolic i/o.

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   The first variable in p is considered as the parameter in the
  --   series of the coefficients in the polynomial on return.
  --   For example, t^3*x + 2*t*x becomes (2*t + t^3)*x.
  --   The number of variables in the polynomial on return is one less
  --   the number of variables in the input polynomial p.
  --   Extra output is written to screen if verbose is true.

  function Series_Polynomial_to_Polynomial
             ( s : Standard_Series_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Converts a polynomial s with coefficients as series
  --   into a polynomial with complex coefficients.
  --   The first variable in the polynomial on return is
  --   the parameter in the series coefficient of s.
  --   Extra output is written to screen if verbose is true.

end Series_and_Polynomials;
