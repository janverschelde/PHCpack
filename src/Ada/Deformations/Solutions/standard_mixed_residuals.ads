with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;

package Standard_Mixed_Residuals is

-- DESCRIPTION :
--   A mixed residual of a polynomial system mixes relative and absolute
--   criteria to decide whether a point is a solution of the system.
--   The computations happen in standard double precision.

  function AbsVal ( z : Vector ) return Vector;

  -- DESCRIPTION :
  --   The i-th component of the vector on return contains the
  --   absolute value (or radius) of z(i).

  function AbsVal ( t : Standard_Complex_Polynomials.Term )
                  return Standard_Complex_Polynomials.Term;
  function AbsVal ( t : Standard_Complex_Laurentials.Term )
                  return Standard_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Returns the same term as t, but with a coefficient
  --   equal to the absolute value of t.cf.

  function AbsVal ( p : Standard_Complex_Polynomials.Poly )
                  return Standard_Complex_Polynomials.Poly;
  function AbsVal ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function AbsVal ( p : Standard_Complex_Laurentials.Poly )
                  return Standard_Complex_Laurentials.Poly;
  function AbsVal ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the polynomials with the same terms as p, but with 
  --   coefficients equal to their absolute values.

  function Residual ( pol,abp : Standard_Complex_Polynomials.Poly;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Laurentials.Poly;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).

  function Residual ( pol,abp : Standard_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_float;
  function Residual ( pol,abp : Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).

end Standard_Mixed_Residuals;
