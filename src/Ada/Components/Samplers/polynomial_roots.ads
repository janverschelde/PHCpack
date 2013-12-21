with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;

package Polynomial_Roots is

-- DESCRIPTION :
--   This package implements a polynomial root finder based on continuation.
--   Below we list all prototypes of the routines to manipulate polynomials
--   in one variable, to set up the homotopy and to track the paths.

  function Eval ( p : Vector; x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the function evaluation of p at x, using Horner's method.

  function Diff ( p : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector of the derivative of p.

  function Roots_of_Unity ( n : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the n roots of unity, i.e.: roots of x^n - 1 = 0.

  procedure Scale ( s0,s1 : in out Complex_Number );

  -- DESCRIPTION :
  --   Applies the scaling s1 = s1/|s1| and s1 = s0/|s1|, if |s1| > 1.

  procedure Scale ( p : in out Vector; s0 : in Complex_Number );

  -- DESCRIPTION :
  --   The i-th coefficient of p is multiplied with s0^(d-i),
  --   where d = degree(p).

  function Scale ( p : Vector; s0 : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   The i-th coefficient of p is multiplied with s0^(d-i),
  --   where d = degree(p).

  procedure Newton ( p,dp : in Vector; x : in out Complex_Number;
                     eps : in double_float; maxit : in natural32;
                     nit : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to the polynomial p, with its derivative
  --   in dp, starting at the point x, until |dx| <= eps or |p(x)| < eps,
  --   not allowing more than maxit iterations.

  function Homotopy ( p,q : Vector; a,t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   Returns (1-t)*p + a*t*q, as t -> 0, the homotopy goes to p.

  procedure Affine_Path_Tracker
              ( p,q : in Vector; a : in Complex_Number;
                s : in out Complex_Number; nbsteps : out natural32 );

  -- DESCRIPTION :
  --   Tracks one path in affine space, using p,q and a as in Homotopy above,
  --   starting at the root s of q.

  procedure Projective_Path_Tracker
              ( p,q : in Vector; a : in Complex_Number;
                s0,s1 : in out Complex_Number; nbsteps : out natural32 );

  -- DESCRIPTION :
  --   Tracks on path in projective space, using p,q and a as in Homotopy 
  --   above, starting at the roots s1/s0 of q.

  procedure Affine_Continuation
              ( p,q : in Vector; a : in Complex_Number; s : in out Vector );
  procedure Affine_Continuation
              ( file : in file_type; p,q : in Vector; a : in Complex_Number;
                s : in out Vector );

  -- DESCRIPTION :
  --   Tracks the paths starting at the roots in s of q, using p,q, and a
  --   as in Homotopy above, in affine coordinates.
  
  procedure Projective_Continuation
              ( p,q : in Vector; a : in Complex_Number;
                s0,s1 : in out Vector );
  procedure Projective_Continuation
              ( file : in file_type; p,q : in Vector; a : in Complex_Number;
                s0,s1 : in out Vector );

  -- DESCRIPTION :
  --   Tracks the paths starting at the roots of s of q, using p,q, and a
  --   as in Homotopy above, in projective coordinates.

  procedure Affine_Solve ( p : in Vector; s : out Vector;
                           maxres : out double_float );
  procedure Affine_Solve ( file : in file_type; p : in Vector;
                           s : out Vector; maxres : out double_float );

  -- DESCRIPTION :
  --   Computes the roots of p in affine coordinates.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   p         coefficient vector of the polynomial.

  -- ON RETURN :
  --   s         the roots of p;
  --   maxres    maximal residual at the roots.

  procedure Projective_Solve
               ( p : in Vector; s0,s1 : out Vector;
                 maxres : out double_float );
  procedure Projective_Solve
               ( file : in file_type; p : in Vector;
                 s0,s1 : out Vector; maxres : out double_float );

  -- DESCRIPTION :
  --   Computes the roots of p in projective coordinates.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   p         coefficient vector of the polynomial.

  -- ON RETURN :
  --   s0,s1     the affine roots of p are in s1(i)/s0(i);
  --   maxres    maximal residual at the roots.

  function Maximal_Residual ( p,s : Vector ) return double_float;
  function Maximal_Residual ( p,s0,s1 : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximal value of the residual, either in affine
  --   or projective coordinates.

  procedure Test_Affine_Roots ( file : in file_type; p,s : in Vector;
                                maxres : out double_float );
  procedure Test_Projective_Roots ( file : in file_type; p,s0,s1 : in Vector;
                                    maxres : out double_float );

  -- DESCRIPTION :
  --   Writes the roots and the residuals on file,
  --   returns in maxres the value of the maximal residual.

end Polynomial_Roots;
