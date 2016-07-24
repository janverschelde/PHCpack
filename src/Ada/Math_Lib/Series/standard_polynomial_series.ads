with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Dense_Series;
with Standard_Series_Polynomials;

package Standard_Polynomial_Series is

-- DESCRIPTION :
--   A series polynomial is a polynomial which has series for coefficients.
--   A polynomial series is a series which has polynomials for coefficients.

  type Poly is record
    deg : integer32; -- the highest power in the series
    cff : Standard_Complex_Poly_Systems.Poly_Sys
            (0..Standard_Dense_Series.max_deg);
     -- coefficients of the series are polynomials in several variables
     -- the number of variables should be the same for all polynomials
  end record;

-- CONSTRUCTORS :

  function Create ( p : Standard_Series_Polynomials.Poly )
                  return Standard_Polynomial_Series.Poly;

  -- DESCRIPTION :
  --   Converts a polynomial p with series coefficients into a series
  --   which has polynomials for coefficients.

  function Create ( p : Standard_Polynomial_Series.Poly )
                  return Standard_Series_Polynomials.Poly;

  -- DESCRIPTION :
  --   Converts a series p which has polynomials for coefficients
  --   into a polynomial with series coefficients.

-- DESTRUCTOR :

  procedure Clear ( p : in out Standard_Polynomial_Series.Poly );

  -- DESCRIPTION :
  --   Deallocates all coefficients in the series p.
  --   On return, p.deg = -1.

end Standard_Polynomial_Series;
