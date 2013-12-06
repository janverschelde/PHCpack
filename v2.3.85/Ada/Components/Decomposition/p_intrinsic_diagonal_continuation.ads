with text_io;                           use text_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package P_Intrinsic_Diagonal_Continuation is

-- DESCRIPTION :
--   While the intrinsic sampling machine suffices to follow
--   paths defined by a diagonal homotopy, efficient evaluations
--   should avoid the explicit doubling of variables.

-- OLD ROUTINES :

  procedure Silent_Diagonal_Continuation
              ( n : in natural; p,q : in Poly;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List );
  procedure Reporting_Diagonal_Continuation
              ( file : in file_type; n : in natural; p,q : in Poly;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Follows the paths moving the affine plane from start to target,
  --   to find witness points on the intersection of p and q.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics when reporting;
  --   n        number of independent variables before doubling;
  --   p,q      polynomials in n variables;
  --   start_b  offset vector of start plane, of range 1..2*n;
  --   target_b offset vector of target plane, of range 1..2*n;
  --   start_v  directions of start plane, of range 1..2;
  --   target_v directions of target plane, of range 1..2;
  --   sols     start solutions in intrinsic coordinates, represented
  --            by coefficients of the directions in start_v.

  -- ON RETURN :
  --   sols     witness points on the intersection of p and q,
  --            on the affine plane defined by target_b and target_v.

-- MORE EFFICIENT ROUTINES :

  function Create ( n : natural; p,q : Poly ) return Eval_Jaco_Mat;

  -- DESCRIPTION :
  --   Returns the evaluable form of the Jacobian matrix of the two
  --   polynomials p and q in n variables.

  procedure Silent_Diagonal_Continuation
              ( n : in natural; p,q : in Eval_Poly; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List );
  procedure Reporting_Diagonal_Continuation
              ( file : in file_type;
                n : in natural; p,q : in Eval_Poly; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Follows the paths moving the affine plane from start to target,
  --   to find witness points on the intersection of p and q.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics when reporting;
  --   n        number of independent variables before doubling;
  --   p,q      polynomials in n variables;
  --   jm       Jacobian Matrix of the system of two polynomials p, q
  --   start_b  offset vector of start plane, of range 1..2*n;
  --   target_b offset vector of target plane, of range 1..2*n;
  --   start_v  directions of start plane, of range 1..2;
  --   target_v directions of target plane, of range 1..2;
  --   sols     start solutions in intrinsic coordinates, represented
  --            by coefficients of the directions in start_v.

  -- ON RETURN :
  --   sols     witness points on the intersection of p and q,
  --            on the affine plane defined by target_b and target_v.

end P_Intrinsic_Diagonal_Continuation;
