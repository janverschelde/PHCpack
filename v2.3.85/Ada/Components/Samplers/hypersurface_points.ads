with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;

package Hypersurface_Points is

-- DESCRIPTION :
--   This package offers facilities to trace point on a hypersurface.

  function Affine_Eval ( p : Eval_Poly; b,v : Vector; t : Complex_Number )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the function p at b + t*v.

  function Projective_Eval ( p : Eval_Poly; b,v : Vector; t : Complex_Number;
                             z : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the multi-homogeneneous function p at b + t*v,
  --   with in z the values for the scaling variables.

  procedure Scale ( b,v : in Vector; t : in Complex_Number;
                    x : out Vector; z : in out Vector );

  -- DESCRIPTION :
  --   Determines the values of the scaling variables z such that
  --   every component of x = b + t*v has modulus one (or less).

  procedure Multihomogenization_Symbols ( n : natural32 );

  -- DESCRIPTION :
  --   Creates a symbol table with x1,..,xn as affine and z1,..,zn
  --   as the extra scaling variables in the multihomogenization.

  function Multihomogenize ( n : natural32; p : Poly ) return Poly;

  -- DESCRIPTION :
  --   Returns the n-homogeneous polynomial, where n is the number of
  --   variables of p.

  procedure Affine_Newton_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  ft,dt : out Complex_Number );
  procedure Projective_Newton_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  z : in Vector; ft,dt : out Complex_Number );

  -- DESCRIPTION
  --   Does one Newton step on the line b + t*v intersected with
  --   the hypersurface p(0) to determine a new value for t.
  --   Applies to moving the line b + t*v on a fixed surface.

  -- ON ENTRY :
  --   p          Horner forms of multivariate polynomials :
  --                p(0) is the original polynomial,
  --                p(i) is the i-th derivative;
  --   b          base point for the line;
  --   v          direction of the line x(t) = b + t*v;
  --   t          parameter to determine location on the line;
  --   z          fixed values of homogenization parameters.

  -- ON RETURN :
  --   t          new location on the line, closer to surface;
  --   ft         previous function value at t;
  --   dt         increment for t, measures also closeness to surface.

  procedure Affine_Newton_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  ft,dt : out Complex_Number );
  procedure Projective_Newton_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  z : in Vector; ft,dt : out Complex_Number );

  -- DESCRIPTION :
  --   Does one Newton step on the line b + t*v intersected with
  --   the moving hypersurface (1-s)*p(0) + s*q(0).

  -- ON ENTRY :
  --   p,q        Horner forms of multivariate polynomials :
  --                p(0) and q(0) are the original polynomials,
  --                p(i) and q(i) are the i-th derivatives;
  --   s          value for the continuation parameter;
  --   b          base point for the line;
  --   v          direction of the line x(t) = b + t*v;
  --   t          parameter to determine location on the line;
  --   z          fixed values of homogenization parameters.

  -- ON RETURN :
  --   t          new location on the line, closer to surface;
  --   ft         previous function value at t;
  --   dt         increment for t, measures also closeness to surface.

  procedure Affine_Correct_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Projective_Correct_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Affine_Correct_for_Line
                ( file : in file_type; p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Projective_Correct_for_Line
                ( file : in file_type; p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to move solution closer to fixed surface.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   p          Horner forms of multivariate polynomials :
  --                p(0) is the original polynomial,
  --                p(i) is the i-th derivative;
  --   b          base point for the line;
  --   v          direction of the line x(t) = b + t*v;
  --   t          parameter to determine location on the line;
  --   z          values for the scaling parameters;
  --   maxit      maximal number of Newton iterations;
  --   eps        desired accuracy.

  -- ON RETURN :
  --   t          new location on the line, closer to surface;
  --   z          updated values for the scaling parameters;
  --   numit      number of iterations performed;
  --   fail       true if desired accuracy not reached.

  procedure Affine_Correct_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Projective_Correct_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Affine_Correct_for_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );
  procedure Projective_Correct_for_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to move solution closer to moving surface.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   p,q        Horner forms of multivariate polynomials :
  --                p(0) and q(0) are the original polynomials,
  --                p(i) and q(i) are the i-th derivatives;
  --   s          continuation parameter in (1-s)*p + s*q, s : 1 -> 0;
  --   b          base point for the line;
  --   v          direction of the line x(t) = b + t*v;
  --   t          parameter to determine location on the line;
  --   z          values for the scaling parameters;
  --   maxit      maximal number of Newton iterations;
  --   eps        desired accuracy.

  -- ON RETURN :
  --   t          new location on the line, closer to surface;
  --   numit      number of iterations performed;
  --   fail       true if desired accuracy not reached.

  procedure Affine_Start_Hypersurface
                ( n,d,k : in natural32; q : out Poly; b,v,t : out Vector );

  -- DESCRIPTION :
  --   Returns a hypersurface of degree d in n variables with a special
  --   random line so that the d roots are roots of unity.

  -- ON ENTRY :
  --   n          number of variables;
  --   d          total degree of the hypersurface;
  --   k          index of the principal variable.

  -- ON RETURN :
  --   q          start polynomial x(k)^d - 1;
  --   b          basis vector with b(k) = 0;
  --   v          direction with v(k) = 0, i /= k;
  --   t          roots of unity.

  procedure Degree_Start_Hypersurface
                ( deg : in Standard_Natural_Vectors.Vector; d : in natural32;
                  v : in Vector; q : out Poly; t : out Vector );

  -- DESCRIPTION :
  --   Returns the start polynomial q(x) = x^deg - 1 = 0 with d solutions
  --   in t to q(v*t) = 0 for the given vector v.

  procedure Affine_Track_Moving_Line
                ( p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean );
  procedure Projective_Track_Moving_Line
                ( p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean );
  procedure Affine_Track_Moving_Line
                ( file : in file_type;
                  p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean );
  procedure Projective_Track_Moving_Line
                ( file : in file_type;
                  p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Performs path tracking from b0 + t*v0 to b1 + t*v1.

  -- ON ENTRY :
  --   file       for diagnostics and intermediate output;
  --   p          polynomial with its derivatives;
  --   (b0,v0)    start line b0 + t*v0;
  --   (b1,v1)    target line b1 + t*v1;
  --   t          value for t at the start of the path;
  --   z          values for scalers in multi-homogenization.
 
  -- ON RETURN :
  --   t          value for t at the end of the path;
  --   z          scaling values used in the homogenization.
  --   numbsteps  number of predictor-corrector steps executed;
  --   fail       true if end of path not reached.

  procedure Affine_Track_Moving_Surface
                ( p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean );
  procedure Projective_Track_Moving_Surface
                ( p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number; z : in out Vector; 
                  numbsteps : out natural32; fail : out boolean );
  procedure Affine_Track_Moving_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean );
  procedure Projective_Track_Moving_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Performs path tracking from the surface q to p.

  -- ON ENTRY :
  --   file       for diagnostics and intermediate output;
  --   p          target polynomial with its derivatives;
  --   q          start polynomial with its derivatives;
  --   (b,v)      fixed line b + t*v;
  --   t          value for t at the start of the path.
 
  -- ON RETURN :
  --   t          value for t at the end of the path;
  --   numbsteps  number of predictor-corrector steps executed;
  --   fail       true if end of path not reached.

  function Maximal_Affine_Residual
                ( p : Eval_Poly; b,v,roots : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximal residual of p(b+roots(i)*v), with one line
  --   output for every root: coordinates and residual.

  function Maximal_Projective_Residual
                ( p : Eval_Poly; b,v,roots : Vector; z : VecVec )
                return double_float;

  -- DESCRIPTION :
  --   Returns the maximal residual of p(b+roots(i)*v) in projective
  --   coordinates.  For every root, this function writes the t-coordinate,
  --   the values for the scalers and the residual.

end Hypersurface_Points;
