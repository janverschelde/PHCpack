with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;

package Standard_Truncated_Series is

-- DESCRIPTION :
--   This packages provides operations to manipulate a truncated
--   power series at a regular point, with standard complex coefficients.
--   The regularity implies that the exponents are natural numbers,
--   and the series is of the form c(0) + c(1)*t + c(2)*t^2 + O(t^3).
--   Thus we identify the truncation with c(0), c(1), c(2),
--   stored as a complex vector of fixed size, with start index at zero.
--   Because of this identification, the addition (and subtraction)
--   of two series are defined by the "+" and "-" of the Vectors package.

  function "*" ( a,b : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the truncated series c = a*b.
 
  -- REQUIRED : a'last = b'last.

  function Inverse ( c : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the inverse of the series, provided c(0) /= 0.

  function "/" ( a,b : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the truncated series c = a/b.

  function Eval ( c : Vector; t : double_float ) return Complex_Number;
  function Eval ( c : Vector; t : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value c[0] + c[1]*t + .. + c[c'last]*t^c'last.

  function sqrt ( c : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the square root of the truncated series with coefficients
  --   in c.  If x = sqrt(c), then x*x = c.

end Standard_Truncated_Series;
