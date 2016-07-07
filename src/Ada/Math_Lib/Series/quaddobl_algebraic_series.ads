with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with QuadDobl_Dense_Series;              use QuadDobl_Dense_Series;

package QuadDobl_Algebraic_Series is

-- DESCRIPTION :
--   An algebraic number is a root of a polynomial in one variable.
--   Similarly, we define an algebraic power series as the root of
--   a polynomial in one variable with coefficients as power series.
--   The calculations are performed in quad double precision.

  function sqrt ( c : Series; i : natural32;
                  verbose : boolean := false ) return Series;

  -- DESCRIPTION :
  --   Applies Newton's method to x^2 - c = 0,
  --   starting at the i-th square root of the zero-th degree.
  --   The degree of the series on return equals c.deg.
  --   If verbose, then the Newton updates dx are written to screen.

  function Root ( c : Series; n,i : natural32;
                  verbose : boolean := false ) return Series;

  -- DESCRIPTION :
  --   Applies Newton's method to x^n - c = 0,
  --   starting at the i-th square root of the zero-th degree.
  --   The degree of the series on return equals c.deg.
  --   If verbose, then the Newton updates dx are written to screen.

end QuadDobl_Algebraic_Series;
