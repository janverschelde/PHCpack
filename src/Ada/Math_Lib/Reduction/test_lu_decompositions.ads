with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with Triple_Double_Numbers;            use Triple_Double_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with Penta_Double_Numbers;             use Penta_Double_Numbers;
with Octo_Double_Numbers;              use Octo_Double_Numbers;
with Deca_Double_Numbers;              use Deca_Double_Numbers;
with Hexa_Double_Numbers;              use Hexa_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Matrices;
with Double_Double_Matrices;
with Triple_Double_Matrices;
with Quad_Double_Matrices;
with Penta_Double_Matrices;
with Octo_Double_Matrices;
with Deca_Double_Matrices;
with Hexa_Double_Matrices;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with TripDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with PentDobl_Complex_Matrices;
with OctoDobl_Complex_Matrices;
with DecaDobl_Complex_Matrices;
with HexaDobl_Complex_Matrices;

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
                A,LU : in Triple_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in triple_double; output : in boolean;
                maxerr : out triple_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Quad_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Penta_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in penta_double; output : in boolean;
                maxerr : out penta_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Octo_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in octo_double; output : in boolean;
                maxerr : out octo_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Deca_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in deca_double; output : in boolean;
                maxerr : out deca_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Hexa_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in hexa_double; output : in boolean;
                maxerr : out hexa_double; fail : out boolean );
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
                A,LU : in TripDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in triple_double; output : in boolean;
                maxerr : out triple_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in QuadDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in PentDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in penta_double; output : in boolean;
                maxerr : out penta_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in OctoDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in octo_double; output : in boolean;
                maxerr : out octo_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in DecaDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in deca_double; output : in boolean;
                maxerr : out deca_double; fail : out boolean );
  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in HexaDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in hexa_double; output : in boolean;
                maxerr : out hexa_double; fail : out boolean );

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
