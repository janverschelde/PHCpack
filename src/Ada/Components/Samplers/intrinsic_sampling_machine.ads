with text_io;                           use text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Intrinsic_Sampling_Machine is

-- DESCRIPTION :
--   This package implements the sampling of a positive dimensional
--   solution component using a parametric description of the linear
--   slicing planes, this is the so-called intrinsic format.

  function Eval ( p : Eval_Poly_Sys; v : VecVec; c : Vector ) return Vector;
  function Eval ( p : Eval_Poly_Sys; b : Vector; v : VecVec; c : Vector )
                return Vector;

  -- DESCRIPTION :
  --   Evaluates the polynomial system at the point x, where x is the
  --   sum of b and c(i)*v(i), for i in v'range.

  -- REQUIRED : v'length = c'length, b'length = v(i)'length, and
  --   b'length is the number of independent variables in p.

  -- ON ENTRY :
  --   p         evaluable form of a polynomial system;
  --   v         directions of the k-plane, k is the number of vectors in v;
  --   b         offset vector of the affine k-plane,
  --             if missing, then we evaluate in projective coordinates;
  --   c         coefficients of the linear combination of vectors in v.

  -- ON RETURN :
  --   p(b+c*v) is the function value of p at the point on the plane.

  function Diff ( jm : Eval_Jaco_Mat; v : VecVec; c : Vector ) return Matrix;
  function Diff ( jm : Eval_Jaco_Mat; b : Vector; v : VecVec; c : Vector )
                return Matrix;

  -- DESCRIPTION :
  --   Returns the evaluation of the Jacobian matrix at the point x,
  --   where x is the sum of b and c(i)*v(i), for i in v'range.
  --   Note that the Jacobian matrix on return contains the values of
  --   the derivatives with respect to c, obtained by the chain rule.

  procedure Moving_Plane
                ( start_b,target_b : in Vector; start_v,target_v : in VecVec;
                  t : in Complex_Number; b : out Vector; v : in out VecVec );

  -- DESCRIPTION :
  --   A family of moving planes can be realized by convex combinations
  --   of the basis and direction vectors of a start and target plane.
  --   This routine computes one plane in this family of moving planes,
  --   given the factor t in the linear pencil of planes.

  -- ON ENTRY :
  --   start_b    is the basis vector of the affine start plane;
  --   target_b   is the basis vector of the affine target plane;
  --   start_v    is the vector of directions of the affine start plane;
  --   target_v   is the vector of directions of the affine target plane;
  --   t          continuation parameter, when t = 0, we obtain the start
  --              plane, at t = 1, we arrive at the target configuration.

  -- ON RETURN :
  --   b          basis point of moving plane : (1-t)*start_b + t * target_b;
  --   v          directions of moving plane : (1-t)*start_v + t*target_v.

  -- REQUIRED :
  --   Memory has already been allocated for the output vector v.
  --   Moreover, v'range = start_v'range = target_v'range; and
  --             b'range = start_b'range = target_b'range.

  procedure Affine_Newton_Refinement
                ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                  b : in Vector; v : in VecVec; c : in out Vector;
                  epsxa,epsfa : in double_float; maxit : in natural );

  -- DESCRIPTION :
  --   Implementation of Newton's method for a system defined on an
  --   affine plane given in its parametric representation by b and v.

  procedure Affine_Solution_Evaluation
                ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                  sol,b : Vector; v : in VecVec );

  -- DESCRIPTION :
  --   The solution vector is decomposed in the parametric representation
  --   of the affine plane and evaluated in the system.

  procedure Silent_Affine_Sampler
                ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                  start_b,target_b : in Vector; start_v,target_v : in VecVec;
                  sols : in out Solution_List );
  procedure Reporting_Affine_Sampler
                ( file : in file_type;
                  p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                  start_b,target_b : in Vector; start_v,target_v : in VecVec;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   Sampling with or without intermediate output between two affine 
  --   planes, from start to target.

  -- ON ENTRY :
  --   file       for intermediate output and results; 
  --   eps        evaluable form of system in original coordinates;
  --   ejm        Jacobi matrix of the system in its original coordinates;
  --   start_b    basis point of the start affine plane;
  --   target_b   basis point of the target affine plane;
  --   start_v    direction vectors of the start affine plane;
  --   target_v   direction vectors of the target affine plane;
  --   sols       contains start solutions as the coefficients of the linear
  --              combination of the directions in start_v.

  -- ON RETURN :
  --   sols       new solutions as coefficient of the linear combination
  --              of the directions in target_v.

end Intrinsic_Sampling_Machine;
