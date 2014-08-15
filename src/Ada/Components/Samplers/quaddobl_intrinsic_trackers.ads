with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Poly_SysFun;      use QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;    use QuadDobl_Complex_Jaco_Matrices;
with Continuation_Parameters;           use Continuation_Parameters;
with QuadDobl_Continuation_Data;        use QuadDobl_Continuation_Data;

package QuadDobl_Intrinsic_Trackers is

-- DESCRIPTION :
--   An intrinsic tracker is a method to follow one path on a solution
--   component of a polynomial system in intrinsic coordinates,
--   in double double precision.
--   A tracker traces a solution along a path from one linear space
--   to another, from the start A to the target B.

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_Affine_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_Projective_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info; k : in natural32;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_Affine_LU_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_Projective_LU_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info; k : in natural32;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_QR_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_QR_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  -- DESCRIPTION :
  --   Tracks the solution from the start to the target plane,
  --   using LU or QR to solve the linear systems.
  --   LU assumes complete intersection, while QR is more general.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   f         a polynomial system in evaluable form, must be
  --             in homogeneous coordinates for projective versions;
  --   jf        Jacobi matrix of the polynomial system f;
  --   s         intrinsic coordinates of solution at start A;
  --   k         scaling index to the largest solution component;
  --   pp        parameters for predictor;
  --   cp        parameters for corrector.

  -- ON RETURN :
  --   s         intrinsic coordinates of solution at target B.

-- GENERIC VERSIONS with POLYNOMIAL FUNCTIONS :

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Silent_LU_Track
               ( ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Reporting_LU_Track
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Silent_QR_Track
               ( ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Reporting_QR_Track
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  -- DESCRIPTION :
  --   Tracks the solution from the start to the target plane,
  --   using LU or QR to solve the linear systems.
  --   LU assumes complete intersection, while QR is more general.
  --   Instead of the evaluable forms of the system f or
  --   the Jacobi matrix jf, these routines take functions.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   ne        number of equations in the polynomial system f;
  --   nv        number of variables in the polynomial system f;
  --   s         intrinsic coordinates of solution at start A;
  --   pp        parameters for predictor;
  --   cp        parameters for corrector.

  -- ON RETURN :
  --   s         intrinsic coordinates of solution at target B.

-- VERSIONS USING LOCAL COORDINATES :

  procedure Silent_Local_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 p : in Matrix; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean );
  procedure Reporting_Local_LU_Track
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 p : in Matrix; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Tracks one solution using local intrinsic coordinates with LU,
  --   moving in the direction to minimize the distance to the target p.

  -- REQUIRED :
  --   p'range(1) = 1..n, n = ambient dimension,
  --   p'range(2) = 0..k, k = codimension of solution set,
  --   p(:,0) is offset vector, vectors in columns 1..k are orthonormal.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   f         a polynomial system made square;
  --   jf        corresponding Jacobian matrix of f;
  --   p         a target k-plane, as offset and orthonormal basis;
  --   s         extrinsic coordinates of a solution;
  --   pp        parameters for predictor;
  --   cp        parameters for corrector.

  -- ON RETURN :
  --   s         extrinsic coordinates of solution at target plane p;
  --   fail      true if target not reached within accuracy requirements,
  --             false if target reached accurately enough.

  procedure Silent_Recentered_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info; pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean );
  procedure Reporting_Recentered_LU_Track
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info; pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Tracks one solution using local intrinsic coordinates with LU,
  --   from start to target plane.

  -- REQUIRED :
  --   p'range(1) = 1..n, n = ambient dimension,
  --   p'range(2) = 0..k, k = codimension of solution set,
  --   p(:,0) is offset vector, vectors in columns 1..k are orthonormal.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   f         a polynomial system made square;
  --   jf        corresponding Jacobian matrix of f;
  --   start     a start k-plane, as offset and orthonormal basis;
  --   target    a target k-plane, as offset and orthonormal basis;
  --   reoriented is true if only start and target differ only in the last
  --             direction vector, otherwise is false;
  --   s         extrinsic coordinates of a solution;
  --   pp        parameters for predictor;
  --   cp        parameters for corrector.

  -- ON RETURN :
  --   s         extrinsic coordinates of solution at target plane p;
  --   fail      true if target not reached within accuracy requirements,
  --             false if target reached accurately enough.

end QuadDobl_Intrinsic_Trackers;
