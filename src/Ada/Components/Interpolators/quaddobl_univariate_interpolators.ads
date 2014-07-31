with QuadDobl_Complex_Numbers;         use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;         use QuadDobl_Complex_Vectors;

package QuadDobl_Univariate_Interpolators is

-- DESCRIPTION :
--   This package implements Newton interpolation with divided differences
--   for functions in one variable over the complex numbers 
--   with double double arithmetic.

-- CONSTRUCTORS :

  function Create ( x,y : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the vector of divided differences of y'range.

  -- REQUIRED : x'range = y'range = 0..d,
  --   with d = degree of interpolating polynomial.

  -- ON ENTRY :
  --   x        interpolation points;
  --   y        function values, y(k) = f(x(k)).

  -- ON RETURN :
  --   vector of divided differences f, which contains
  --     f(0) = f[x(0)],
  --     f(1) = f[x(0)x(1)], ...
  --     f(k) = f[x(0)x(1)..x(k)], k = 0..d = x'last.
  --   The interpolating polynomial is then represented as the
  --   sum of terms of f(k)*(x-x(1))*(x-x(2))*..*(x-x(k-1)), 
  --   for k ranging from 0 to d.

  function Expand ( f,x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given the divided differences f and the interpolation points x
  --   (same x that was used to create f), Expand returns the coefficients
  --   of the polynomial in expanded form.  If c = Expand(f,x), then c
  --   represents c(0) + c(1)*x + c(2)*x^2 + .. + c(c'last)*x^c'last,
  --   where c'last = x'last.

-- EVALUATORS :

  function Evalf ( f,x : Vector; a : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Given the representation of the interpolating polynomial p by its
  --   divided differences f and the interpolation points x, Evalf returns
  --   the function value of p at the point a.

  function Evalc ( c : Vector; x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Given the representation of the interpolating polynomial p by its
  --   coefficients (created with Expand), the value of p at the point x
  --   is returned.

end QuadDobl_Univariate_Interpolators;
