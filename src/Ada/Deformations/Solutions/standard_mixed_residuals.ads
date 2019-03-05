with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package Standard_Mixed_Residuals is

-- DESCRIPTION :
--   A mixed residual of a polynomial system mixes relative and absolute
--   criteria to decide whether a point is a solution of the system.
--   The computations happen in standard double precision.

  function AbsVal ( z : Vector ) return Vector;

  -- DESCRIPTION :
  --   The i-th component of the vector on return contains the
  --   absolute value (or radius) of z(i).

  function AbsVal ( t : Term ) return Term;

  -- DESCRIPTION :
  --   Returns the same term as t, but with a coefficient
  --   equal to the absolute value of t.cf.

  function AbsVal ( p : Poly ) return Poly;

  -- DESCRIPTION :
  --   Returns the polynomials with the same terms as p, but with 
  --   coefficients equal to their absolute values.

  function Residual ( pol,abp : Poly; z : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial pol at z,
  --   where abp = AbsVal(pol).

end Standard_Mixed_Residuals;
