with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Vectors;
with Triple_Double_Matrices;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;

package Test_TripDobl_QRLS_Solvers is

-- DESCRIPTION :
--   Tests the QR decomposition and Least Squares linear system solving
--   in triple double precision.

  function Extract_Upper_Triangular
                ( a : Triple_Double_Matrices.Matrix )
                return Triple_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

  function Extract_Upper_Triangular
                ( a : TripDobl_Complex_Matrices.Matrix )
                return TripDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

  function Differences ( a,b : in Triple_Double_Matrices.Matrix )
                       return triple_double;

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

  function Differences ( a,b : in TripDobl_Complex_Matrices.Matrix )
                       return triple_double;

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

  function Orthogonality_Check_Sum
             ( q : Triple_Double_Matrices.Matrix ) return triple_double;

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q 
  --   with all its following columns.

  function Orthogonality_Check_Sum 
             ( q : TripDobl_Complex_Matrices.Matrix ) return triple_double;

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q with
  --   its following columns.

  procedure Test_QRD ( a,q,r : in Triple_Double_Matrices.Matrix;
                       output : in boolean );

  -- DESCRIPTION :
  --   Given in q and r is the QR decomposition of the real matrix a.
  --   Tests the properties of q and r with output if output.

  procedure Test_QRD ( a,q,r : in TripDobl_Complex_Matrices.Matrix;
                       output : in boolean );

  -- DESCRIPTION :
  --   Given in q and r is the QR decomposition of the complex matrix a.
  --   Tests the properties of q and r with output if output.

  procedure TripDobl_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Triple_Double_Matrices.Matrix;
                b : in Triple_Double_Vectors.Vector;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a real n-by-m matrix a and an n-vector b,
  --   tests the least squares solver with pivoting if piv
  --   and with output of all numbers if output.

  procedure TripDobl_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Triple_Double_Matrices.Matrix;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a real n-by-m matrix,
  --   tests the QR decomposition with pivoting if piv
  --   and with output of all numbers if output.

  procedure TripDobl_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for an n-by-m real matrix and runs a test
  --   on the QR decomposition with pivoting if piv is true.

  procedure TripDobl_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for an n-by-m real matrix and runs a test
  --   on the least squares solver with pivoting if piv is true.

  procedure TripDobl_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for the number of tests on the QR decomposition
  --   of random n-by-m real matrices with pivoting if piv is true.

  procedure TripDobl_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for the number of tests on the least squares solver
  --   of random n-by-m real matrices with pivoting if piv is true.

  procedure TripDobl_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : TripDobl_Complex_Matrices.Matrix;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a complex n-by-m matrix,
  --   tests the QR decomposition with pivoting if piv
  --   and with output of all numbers if output.

  procedure TripDobl_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : TripDobl_Complex_Matrices.Matrix;
                b : TripDobl_Complex_Vectors.Vector;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a complex n-by-m matrix a and an n-vector b,
  --   tests the least squares solver with pivoting if piv
  --   and with output of all numbers if output.

  procedure TripDobl_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a complex n-by-m matrix and tests the QR decomposition
  --   with pivoting if piv is true.

  procedure TripDobl_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a complex n-by-m matrix and complex n-dimensional
  --   right hand side vector and tests the least squares solver,
  --   with pivoting if piv is true.

  procedure TripDobl_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a number of tests and times the QR decomposition
  --   on n-by-m matrices with pivoting if piv is true.

  procedure TripDobl_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a number of tests and times the least squares
  --   solver on n-by-m matrices with pivoting if piv is true.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_TripDobl_QRLS_Solvers;
