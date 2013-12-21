with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;

package Standard_Intrinsic_Newton is

-- DESCRIPTION :
--   Implementation of Newton's method for solutions of polynomial
--   systems restricted to linear spaces.

  function Affine_Eval
             ( jf : Eval_Jaco_Mat; p : Matrix; x : Vector )
	     return Matrix;

  -- DESCRIPTION :
  --   Returns the value of the Jacobi matrix at the point given in
  --   affine extrinsic coordinates x with respect to the plane p.

  generic
    with function jf ( x : Vector ) return Matrix;
  function Generic_Affine_Eval ( p : Matrix; x : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the value of the Jacobi matrix at the point given in
  --   affine extrinsic coordinates x with respect to the plane p.
  --   The function jf evaluations the Jacobi matrix at the point.

  function Projective_Eval 
             ( jf : Eval_Jaco_Mat; p : Matrix; x : Vector; k : natural32 )
             return Matrix;

  -- DESCRIPTION :
  --   Returns the value of the Jacobi matrix at the point given in
  --   affine extrinsic coordinates x with respect to the plane p.
  --   The intrinsic coordinates are in projective space.
  --   When k = 0, Projective_Eval returns the same matrix as Affine_Eval,
  --   otherwise, the matrix on return will still have the same range,
  --   but the k-th intrinsic coordinate will be considered as fixed. 

  procedure Affine_LU_Newton
              ( f : in Poly_Sys; 
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Affine_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_LU_Newton
              ( f : in Poly_Sys; 
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  procedure Affine_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Affine_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Projective_LU_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Projective_LU_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );

  procedure Affine_QR_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Affine_QR_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_QR_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_QR_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  procedure Affine_SV_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Affine_SV_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Projective_SV_Newton
              ( f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Projective_SV_Newton
              ( file : in file_type; f : in Poly_Sys;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Implementation of Newton's method, with some variations:
  --     LU : uses plain LU factorization with optional estimate
  --          for condition number, needs complete intersection;
  --     QR : uses QR to deal with overconstrained linear systems;
  --     SV : computes singular values, for special occasions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then the procedure writes nothing;
  --   f        a polynomial system;
  --   p        orthonormal representation of a linear space;
  --   x        approximate solution in intrinsic coordinates;
  --   k        index of coordinate which stays fixed;
  --   epsax    accuracy requirement for absolute error on x;
  --   epsrx    accuracy requirement for relative error on x;
  --   epsaf    accuracy requirement for residual f(x);
  --   epsrf    accuracy requirement for f(x)/x;
  --   maxit    maximal number of iterations allowed.

  -- ON RETURN :
  --   x        more accurate solution in intrinsic coordinates;
  --   incax    estimate for absolute error on x (max norm of increment);
  --   incrx    estimate for relative error on x;
  --   resaf    max norm of absolute residual, i.e.: |f(x)|;
  --   resrf    resaf divided by the max norm of x;
  --   nit      number of iterations used;
  --   rco      estimate for inverse condition number (optional when LU);
  --   sv       vector of singular values returned with SV;
  --   nit      number of Newton iterations used;
  --   fail     true if any of the accuracy requirements is not reached
  --            within maxit steps; false otherwise.

  procedure Affine_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Affine_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  procedure Affine_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Affine_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Projective_LU_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );
  procedure Projective_LU_Newton
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );

  procedure Affine_QR_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Affine_QR_Newton 
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_QR_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );
  procedure Projective_QR_Newton 
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  procedure Affine_SV_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Affine_SV_Newton 
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Projective_SV_Newton
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );
  procedure Projective_SV_Newton 
              ( file : in file_type;
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
		p : in Matrix; x : in out Vector; k : in natural32;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );

  -- DESCRIPTION :
  --   These routines are called by the corresponding routines above.
  --   For repeated application of Newton's method on the same system,
  --   calling these routines directly avoids the repeated creation
  --   of Jacobi matrices and Horner forms of the polynomials.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then the procedure writes nothing;
  --   f        evaluable form a polynomial system;
  --   jf       evaluable form of the Jacobi matrix of f;
  --   p        orthonormal representation of a linear space;
  --   x        approximate solution in intrinsic coordinates;
  --   k        index of coordinate which stays fixed;
  --   epsax    accuracy requirement for absolute error on x;
  --   epsrx    accuracy requirement for relative error on x;
  --   epsaf    accuracy requirement for residual f(x);
  --   epsrf    accuracy requirement for f(x)/x;
  --   maxit    maximal number of iterations allowed.

  -- ON RETURN :
  --   x        more accurate solution in intrinsic coordinates;
  --   incax    estimate for absolute error on x (max norm of increment);
  --   incrx    estimate for relative error on x;
  --   resaf    max norm of absolute residual, i.e.: |f(x)|;
  --   resrf    resaf divided by the max norm of x;
  --   nit      number of iterations used;
  --   rco      estimate for inverse condition number (optional when LU);
  --   sv       vector of singular values returned with SV;
  --   nit      number of Newton iterations used;
  --   fail     true if any of the accuracy requirements is not reached
  --            within maxit steps; false otherwise.

-- GENERIC VERSIONS :

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_Affine_LU_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_Affine_LU_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_Affine_LU_RCO_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_Affine_LU_RCO_Newton
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                rco : out double_float; fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_Affine_QR_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_Affine_QR_Newton 
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                fail : out boolean );

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_Affine_SV_Newton
              ( n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );

  generic 
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_Affine_SV_Newton 
              ( file : in file_type;
                n : in natural32; p : in Matrix; x : in out Vector;
                epsax,epsrx,epsaf,epsrf : in double_float;
                incax,incrx,resaf,resrf : out double_float;
                nit : out natural32; maxit : in natural32;
                sv : out Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Instead of evaluable forms, f and its Jacobi matrix evaluation
  --   should be given as generic functions.
  --   The number of equations in f is an additional parameter.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then the procedure writes nothing;
  --   n        number of equations in f;
  --   p        orthonormal representation of a linear space;
  --   x        approximate solution in intrinsic coordinates;
  --   k        index of coordinate which stays fixed;
  --   epsax    accuracy requirement for absolute error on x;
  --   epsrx    accuracy requirement for relative error on x;
  --   epsaf    accuracy requirement for residual f(x);
  --   epsrf    accuracy requirement for f(x)/x;
  --   maxit    maximal number of iterations allowed.

  -- ON RETURN :
  --   x        more accurate solution in intrinsic coordinates;
  --   incax    estimate for absolute error on x (max norm of increment);
  --   incrx    estimate for relative error on x;
  --   resaf    max norm of absolute residual, i.e.: |f(x)|;
  --   resrf    resaf divided by the max norm of x;
  --   nit      number of iterations used;
  --   rco      estimate for inverse condition number (optional when LU);
  --   sv       vector of singular values returned with SV;
  --   nit      number of Newton iterations used;
  --   fail     true if any of the accuracy requirements is not reached
  --            within maxit steps; false otherwise.

end Standard_Intrinsic_Newton;
