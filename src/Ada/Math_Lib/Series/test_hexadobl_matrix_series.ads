with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Complex_Series_Matrices;

package Test_HexaDobl_Matrix_Series is

-- DESCRIPTION :
--   Tests matrices of truncated dense power series
--   in hexa double precision.

  procedure Write ( v : in HexaDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Basic output of a vector of series in hexa double precision.

  procedure Write ( A : in HexaDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPTION :
  --   Basic output of a matrix of standard series.

  procedure Write ( A : in HexaDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Basic output of a hexa double complex matrix.

  procedure HexaDobl_Random_Linear_Solve ( n,degree : in integer32 );

  -- DESCRIPTION :
  --   Generates a random matrix A and right hand side vector b
  --   of dimension n and solves the linear system A*x = b,
  --   in hexa double precision.
  --   The series are all of the given degree.

  function Trunc ( A : HexaDobl_Complex_Series_Matrices.Matrix )
                 return HexaDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Takes the zero-th degree coefficient of every element in A.

  procedure Zero_QRD ( A : HexaDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows the QR decomposition of the zero degree terms of A.

  procedure Solve_Normal_Equations
              ( A : in HexaDobl_Complex_Series_Matrices.Matrix;
                x,b : in HexaDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations,
  --   in hexa double precision.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

  procedure Test_Normality
              ( qr : in HexaDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedure computes the norm of every column in qr,
  --   in hexa double precision.

  procedure Test_Orthogonality
              ( qr : in HexaDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes all products of the vectors with
  --   all other vectors to test whether they are orthogonal.

  procedure Test_Basis
              ( wrk,A : in HexaDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality, in hexa double precision.

  procedure QR_Solve_Least_Squares
              ( A : in HexaDobl_Complex_Series_Matrices.Matrix;
                x,b : in HexaDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies QR decomposition to the matrix A and then solves the
  --   system A*x = b with this QR decomposition.
  --   The computed solution is compared to the constructed solution
  --   and the residual is computed.

  procedure HexaDobl_Random_Least_Squares ( n,m,degree : in integer32 );

  -- DESCRIPTION :
  --   Generates a random n-by-m matrix A with series of the given degree.
  --   A random solution x is generated and then the right hand side
  --   vector b is A*x.  for n > m this is an overdetermined linear
  --   system A*x = b, where the solution x is the generated vector.
  --   Computations are done in hexa double precision.

  procedure Test_Solving;

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degree of the series, solves a random system.

  procedure HexaDobl_Test_Conversion ( n,m,d : in integer32 );

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

  procedure Test_Conversion;

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and the degree of the series.

  procedure Main;

  -- DESCRIPTION :
  --   Displays the menu of test procedures and runs the selected test.

end Test_HexaDobl_Matrix_Series;
