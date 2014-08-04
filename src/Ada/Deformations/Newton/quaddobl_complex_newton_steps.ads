with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;

package QuadDobl_Complex_Newton_Steps is

-- DESCRIPTION :
--   This package provides two generic routines to execute one step
--   of Newton's method, using a singular value decomposition on the
--   Jacobian matrix, all in standard complex floating-point arithmetic.

  generic

    with function f  ( x : Vector ) return Vector;  -- returns function value
    with function jm ( x : Vector ) return Matrix;  -- returns Jacobi matrix

  procedure Silent_Newton_Step
              ( n : in natural32; z : in out Vector; tol : in double_float;
                err,rco,res : out quad_double; rank : out natural32 );

  -- DESCRIPTION :
  --   Applies one step with Newton's method on the system f(x) = 0,
  --   using standard machine precision complex floating-point arithmetic.

  -- ON ENTRY :
  --   n        number of equations, may be larger than z'length;
  --   z        initial approximation for the root of f;
  --   tol      tolerance to decide the numerical rank.

  -- ON RETURN :
  --   z        refined approximation for the root of f;
  --   err      maximum norm of the correction vector;
  --   rco      estimate for the inverse condition number;
  --   rank     numerical rank of Jacobian matrix at root.

  generic

    with function f  ( x : Vector ) return Vector;  -- returns function value
    with function jm ( x : Vector ) return Matrix;  -- returns Jacobi matrix

  procedure Reporting_Newton_Step
              ( file : in file_type;
                n : in natural32; z : in out Vector; tol : in double_float;
                err,rco,res : out quad_double; rank : out natural32 );

  -- DESCRIPTION :
  --   Applies one step with Newton's method on the system f(x) = 0,
  --   using standard machine precision complex floating-point arithmetic.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of equations, may be larger than z'length;
  --   z        initial approximation for the root of f;
  --   tol      tolerance to decide the numerical rank.

  -- ON RETURN :
  --   z        refined approximation for the root of f;
  --   err      maximum norm of the correction vector;
  --   rco      estimate for the inverse condition number;
  --   res      maximum norm of the residual vector;
  --   rank     numerical rank of the Jacobian matrix.

  generic

    with function f  ( x : Vector ) return Vector;  -- returns function value
    with function jm ( x : Vector ) return Matrix;  -- returns Jacobi matrix

  procedure Silent_Newton_Step_with_Singular_Values
              ( n : in natural32; z : in out Vector; tol : in double_float;
                err,rco,res : out quad_double;
                s : out Vector; rank : out natural32 );

  -- DESCRIPTION :
  --   Applies one step with Newton's method on the system f(x) = 0,
  --   using standard machine precision complex floating-point arithmetic.
  --   This version returns the vector of singular values computed by the
  --   SVD algorithm on the Jacobi matrix at z.

  -- ON ENTRY :
  --   n        number of equations, may be larger than z'length;
  --   z        initial approximation for the root of f;
  --   tol      tolerance to decide the numerical rank.

  -- ON RETURN :
  --   z        refined approximation for the root of f;
  --   err      maximum norm of the correction vector;
  --   rco      estimate for the inverse condition number;
  --   res      maximum norm of the residual vector;
  --   s        vector of range min0(n+1,z'last) with singular values,
  --            of the Jacobian matrix at the z on input;
  --   rank     numerical rank of the Jacobian matrix.

  generic

    with function f  ( x : Vector ) return Vector;  -- returns function value
    with function jm ( x : Vector ) return Matrix;  -- returns Jacobi matrix

  procedure Reporting_Newton_Step_with_Singular_Values
              ( file : in file_type;
                n : in natural32; z : in out Vector; tol : in double_float;
                err,rco,res : out quad_double;
                s : out Vector; rank : out natural32 );

  -- DESCRIPTION :
  --   Applies one step with Newton's method on the system f(x) = 0,
  --   using standard machine precision complex floating-point arithmetic.
  --   This version returns the vector of singular values computed by the
  --   SVD algorithm on the Jacobi matrix at z.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of equations, may be larger than z'length;
  --   z        initial approximation for the root of f;
  --   tol      tolerance to decide the numerical rank.

  -- ON RETURN :
  --   z        refined approximation for the root of f;
  --   err      maximum norm of the correction vector;
  --   rco      estimate for the inverse condition number;
  --   res      maximum norm of the residual vector;
  --   s        vector of range min0(n+1,z'last) with singular values,
  --            of the Jacobian matrix at the z on input;
  --   rank     numerical rank of Jacobian matrix.

end QuadDobl_Complex_Newton_Steps;
