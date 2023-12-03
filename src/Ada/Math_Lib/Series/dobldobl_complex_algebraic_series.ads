with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Series;            use DoblDobl_Complex_Series;

package DoblDobl_Complex_Algebraic_Series is

-- DESCRIPTION :
--   An algebraic number is a root of a polynomial in one variable.
--   Similarly, we define an algebraic power series as the root of
--   a polynomial in one variable with coefficients as power series.

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

  function Poly_Eval ( p : Vector; z : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the value of the polynomial with coefficients in p,
  --   evaluated at the series z.

  function Poly_Diff ( p : Vector; z : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the value of the derivative of the polynomial 
  --   with coefficients in p, evaluated at the series z.

  function Poly_Root ( p : Vector; z0 : Complex_Number; c : Series; 
                       verbose : boolean := false ) return Series;

  -- DESCRIPTION :
  --   Returns the series expansion z(t) of p(x) - c(t) = 0,
  --   applying Newton iteration at a root of p.
  --
  -- ON ENTRY :
  --   p        coefficient vector of a polynomial p in x,
  --            p(i) is the coefficient of x^i of p;
  --   z0       leading coefficient of the z(t): p(z0) = z0;
  --   c        right hand side series for the equation,
  --            c.deg is the degree of the series on return;
  --   verbose  if true, then the Newton updates are written.

end DoblDobl_Complex_Algebraic_Series;
