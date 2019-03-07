with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;

package DoblDobl_Mixed_Residuals is

-- DESCRIPTION :
--   A mixed residual of a polynomial system mixes relative and absolute
--   criteria to decide whether a point is a solution of the system.
--   The computations happen in double double precision.

  function AbsVal ( z : Vector ) return Vector;

  -- DESCRIPTION :
  --   The i-th component of the vector on return contains the
  --   absolute value (or radius) of z(i).

  function AbsVal ( t : DoblDobl_Complex_Polynomials.Term )
                  return DoblDobl_Complex_Polynomials.Term;
  function AbsVal ( t : DoblDobl_Complex_Laurentials.Term )
                  return DoblDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Returns the same term as t, but with a coefficient
  --   equal to the absolute value of t.cf.

  function AbsVal ( p : DoblDobl_Complex_Polynomials.Poly )
                  return DoblDobl_Complex_Polynomials.Poly;
  function AbsVal ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function AbsVal ( p : DoblDobl_Complex_Laurentials.Poly )
                  return DoblDobl_Complex_Laurentials.Poly;
  function AbsVal ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the polynomials with the same terms as p, but with 
  --   coefficients equal to their absolute values.

  function Residual ( pol,abp : DoblDobl_Complex_Polynomials.Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laurentials.Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).

  function Residual ( pol,abp : DoblDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_double;
  function Residual ( pol,abp : DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the mixed residual of the polynomial(s) pol at z,
  --   where abp = AbsVal(pol).

end DoblDobl_Mixed_Residuals;
