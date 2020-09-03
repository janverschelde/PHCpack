with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Vectors;
with Deca_Double_Matrices;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;

package Test_DecaDobl_QRLS_Solvers is

-- DESCRIPTION :
--   Tests the QR decomposition and Least Squares linear system solving
--   in deca double precision.

  function Extract_Upper_Triangular
                ( a : Deca_Double_Matrices.Matrix )
                return Deca_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

  function Extract_Upper_Triangular
                ( a : DecaDobl_Complex_Matrices.Matrix )
                return DecaDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

  function Differences ( a,b : in Deca_Double_Matrices.Matrix )
                       return Deca_Double;

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

  function Differences ( a,b : in DecaDobl_Complex_Matrices.Matrix )
                       return Deca_Double;

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

  function Orthogonality_Check_Sum
             ( q : Deca_Double_Matrices.Matrix ) return Deca_Double;

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q 
  --   with all its following columns.

  function Orthogonality_Check_Sum 
             ( q : DecaDobl_Complex_Matrices.Matrix ) return Deca_Double;

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q with
  --   its following columns.

  procedure Test_QRD ( a,q,r : in Deca_Double_Matrices.Matrix;
                       output : in boolean );

  -- DESCRIPTION :
  --   Given in q and r is the QR decomposition of the real matrix a.
  --   Tests the properties of q and r with output if output.

  procedure Test_QRD ( a,q,r : in DecaDobl_Complex_Matrices.Matrix;
                       output : in boolean );

  -- DESCRIPTION :
  --   Given in q and r is the QR decomposition of the complex matrix a.
  --   Tests the properties of q and r with output if output.

  procedure DecaDobl_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Deca_Double_Matrices.Matrix;
                b : in Deca_Double_Vectors.Vector;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a real n-by-m matrix a and an n-vector b,
  --   tests the least squares solver with pivoting if piv
  --   and with output of all numbers if output.

  procedure DecaDobl_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Deca_Double_Matrices.Matrix;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a real n-by-m matrix,
  --   tests the QR decomposition with pivoting if piv
  --   and with output of all numbers if output.

  procedure DecaDobl_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for an n-by-m real matrix and runs a test
  --   on the QR decomposition with pivoting if piv is true.

  procedure DecaDobl_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for an n-by-m real matrix and runs a test
  --   on the least squares solver with pivoting if piv is true.

  procedure DecaDobl_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for the number of tests on the QR decomposition
  --   of random n-by-m real matrices with pivoting if piv is true.

  procedure DecaDobl_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for the number of tests on the least squares solver
  --   of random n-by-m real matrices with pivoting if piv is true.

  procedure DecaDobl_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : DecaDobl_Complex_Matrices.Matrix;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a complex n-by-m matrix,
  --   tests the QR decomposition with pivoting if piv
  --   and with output of all numbers if output.

  procedure DecaDobl_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : DecaDobl_Complex_Matrices.Matrix;
                b : DecaDobl_Complex_Vectors.Vector;
                output : in boolean );

  -- DESCRIPTION :
  --   On input is a complex n-by-m matrix a and an n-vector b,
  --   tests the least squares solver with pivoting if piv
  --   and with output of all numbers if output.

  procedure DecaDobl_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a complex n-by-m matrix and tests the QR decomposition
  --   with pivoting if piv is true.

  procedure DecaDobl_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a complex n-by-m matrix and complex n-dimensional
  --   right hand side vector and tests the least squares solver,
  --   with pivoting if piv is true.

  procedure DecaDobl_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a number of tests and times the QR decomposition
  --   on n-by-m matrices with pivoting if piv is true.

  procedure DecaDobl_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean );

  -- DESCRIPTION :
  --   Prompts for a number of tests and times the least squares
  --   solver on n-by-m matrices with pivoting if piv is true.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_DecaDobl_QRLS_Solvers;
