with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Matrices;
with Double_Double_Matrices;
with Quad_Double_Matrices;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;

package Test_LU_Decompositions is

-- DESCRIPTION :
--   This packages offers routines to test the accuracy of the
--   LU decomposition for various precision levels.

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Standard_Floating_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean );

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Double_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_double; output : in boolean;
                maxerr : out double_double; fail : out boolean );

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Quad_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean );

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Standard_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean );

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in DoblDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_double; output : in boolean;
                maxerr : out double_double; fail : out boolean );

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in QuadDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean );

  -- DESCRIPTION :
  --   Checks whether P*A = L*U, for a given floating point matrix A,
  --   its LU decomposition in LU, and permutation in ipvt.

  -- ON ENTRY :
  --   n        number of rows and columns of matrices A and LU;
  --   A        original n-by-n matrix;
  --   LU       contains upper triangular part of A and the multipliers
  --            below the diagonal used in the factorization;
  --   ipvt     indices of pivots used, defines permutation matrix;
  --   tol      tolerance on difference between the P*A and L*U;
  --   output   if true, then intermediate output is shown.

  -- ON RETURN :
  --   maxerr   largest error componentwise between P*A and L*U;
  --   fail     true if maxerr is larger than the tolerance tol.

end Test_LU_Decompositions;
