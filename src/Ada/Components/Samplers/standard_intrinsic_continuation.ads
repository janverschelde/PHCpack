with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Continuation_Parameters;           use Continuation_Parameters;
with Standard_Continuation_Data;        use Standard_Continuation_Data;

package Standard_Intrinsic_Continuation is

-- DESCRIPTION :
--   Intrinsic continuation methods use intrinsic coordinates of
--   solutions on a solution component of a polynomial system.
--   The continuation starts at solutions on a start plane A
--   and ends at solutions on a target plane B.
--   The path between A and B is defined by a generic function.

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_Affine_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_Projective_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_Affine_LU_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_Projective_LU_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Silent_QR_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function Path ( t : Complex_Number ) return Matrix;
  procedure Reporting_QR_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  -- DESCRIPTION :
  --   Tracks the solution from the start to the target plane,
  --   using LU or QR to solve the linear systems.
  --   LU assumes complete intersection, while QR is more general.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   f         a polynomial system in evaluable form;
  --   jf        Jacobi matrix of the polynomial system f;
  --   s         intrinsic coordinates of solutions at start plane A;
  --   p         parameters for predictor;
  --   c         parameters for corrector.

  -- ON RETURN :
  --   s         intrinsic coordinates of solutions at target plane B.

  procedure LU_Validate
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  procedure LU_Validate
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  procedure SV_Validate
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  procedure SV_Validate
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );

  -- DESCRIPTION :
  --   Applies a couple of Newton iterations to validate whether
  --   the solutions on the plane satisfy the polynomial system,
  --   using LU (for complete intersections) or SVD (in general).

-- GENERIC VERSIONS with POLYNOMIAL FUNCTIONS :

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Silent_LU_Continue
               ( ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Reporting_LU_Continue
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Silent_QR_Continue
               ( ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  generic           -- A = Path(0) is start, B = Path(1) is target plane
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
    with function Path ( t : Complex_Number ) return Matrix;
  procedure G_Reporting_QR_Continue
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars );

  -- DESCRIPTION :
  --   These are the generic versions of the routines above.
  --   Tracks the solution from the start to the target plane,
  --   using LU or QR to solve the linear systems.
  --   LU assumes complete intersection, while QR is more general.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   ne        number of polynomials in the polynomial sysem f;
  --   nv        number of variables in the polynomial sysem f;
  --   s         intrinsic coordinates of solutions at start plane A;
  --   p         parameters for predictor;
  --   c         parameters for corrector.

  -- ON RETURN :
  --   s         intrinsic coordinates of solutions at target plane B.

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_LU_Validate
               ( n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_LU_Validate
               ( file : in file_type; n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Silent_SV_Validate
               ( n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );
  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Reporting_SV_Validate
               ( file : in file_type; n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 );

  -- DESCRIPTION :
  --   Applies a couple of Newton iterations to validate whether
  --   the solutions on the plane satisfy the polynomial system,
  --   using LU (for complete intersections) or SVD (in general).
  --   The extra variable n is the number of equations in f.

-- using LOCAL COORDINATES :

  procedure Silent_Local_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info_Array;
                 pp : in Pred_Pars; cp : in Corr_Pars );
  procedure Reporting_Local_LU_Continue
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info_Array;
                 pp : in Pred_Pars; cp : in Corr_Pars );

  -- DESCRIPTION :
  --   Local coordinates are applied in a recentering algorithm
  --   to track all paths from the current set of solutions to
  --   lie on the target plane provided by the matrix p.

  -- REQUIRED :
  --   p'range(1) = 1..n, n = ambient dimension,
  --   p'range(2) = 0..k, k = codimension of solution set,
  --   p(:,0) is offset vector, vectors in columns 1..k are orthonormal.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics,
  --             without file, no output is produced;
  --   f         a polynomial system made square;
  --   jf        corresponding Jacobian matrix of f;
  --   start     the start k-plane on which the solutions of f lie;
  --   target    a target k-plane, as offset and orthonormal basis;
  --   reoriented indicates whether start and target differ only in
  --             their last direction vector, is false otherwise;
  --   s         extrinsic coordinates of solutions;
  --   pp        parameters for predictor;
  --   cp        parameters for corrector.

  -- ON RETURN :
  --   s         extrinsic coordinates of solutions at target plane p.

end Standard_Intrinsic_Continuation;
