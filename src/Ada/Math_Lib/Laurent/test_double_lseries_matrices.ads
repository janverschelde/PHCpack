with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecVecs;

package Test_Double_Lseries_Matrices is

-- DESCRIPTION :
--   A matrix of Laurent series is a tuple of
--   1) a matrix of leading exponents, and
--   2) a 3-dimensional vector of vector of vectors
--   with the coefficients of the power series.
--   The procedures in this package test the LU factorization
--   of a matrix of Laurent series, in double precision.

  procedure Matrix_Matrix_Product
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Clead : out Standard_Integer_Matrices.Matrix;
                Ccffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns the product of two square matrices A and B in C.

  -- ON ENTRY :
  --   nrows    number of rows of all matrices;
  --   ncols    number of columns of all matrices;
  --   deg      degree of the series in all matrices;
  --   Alead    leading exponents of the series in the matrix A;
  --   Acffs    coefficients of the series in the matrix A;
  --   Blead    leading exponents of the series in the matrix B;
  --   Bcffs    coefficients of the series in the matrix B.

  -- ON RETURN :
  --   Clead    leading exponents of the product of A with B;
  --   Ccffs    coefficients of the product of A with B.

  procedure Write_Difference
              ( deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Writes the difference between the matrices A and B.

  procedure Plain_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Llead : out Standard_Integer_Matrices.Matrix;
                Lcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Ulead : out Standard_Integer_Matrices.Matrix;
                Ucffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Applies the plain LU factorization, plain means without pivoting,
  --   factoring the matrix A as a product of a lower triangular L with
  --   an upper triangular matrix U.  For testing purposes, the matrices
  --   L and U are returned explicitly.

  -- ON ENTRY :
  --   nrows    number of rows of the matrix;
  --   ncols    number of columns of the matrix;
  --   deg      degree of the series in all matrices;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series.

  -- ON RETURN :
  --   Llead    leading exponents of the lower triangular factor L
  --            that has ones on its diagonal;
  --   Lcffs    coefficients of the lower triangular factor L;
  --   Ulead    leading exponents of the upper triangular factor U;
  --   Ucffs    coefficients of the upper triangular factor U.

  procedure Test_Plain_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Tests the LU factorization on a general, random matrix.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series.

  procedure Lower_Triangular_Part
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Llead : out Standard_Integer_Matrices.Matrix;
                Lcffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns the lower triangular part of A, with ones on the diagonal.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series in the matrix.

  -- ON RETURN :
  --   Llead    leading exponents of the lower triangular part of A;
  --   Lcffs    coefficients of the series in the lower triangular part.

  procedure Upper_Triangular_Part
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Ulead : out Standard_Integer_Matrices.Matrix;
                Ucffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns the upper triangular part of A,
  --   which includes the diagonal.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series in the matrix.

  -- ON RETURN :
  --   Ulead    leading exponents of the upper triangular part of A;
  --   Ucffs    coefficients of the series in the upper triangular part.

  procedure Permute
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : in Standard_Integer_Vectors.Vector;
                Blead : out Standard_Integer_Matrices.Matrix;
                Bcffs : out Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns in B the pivoting applied to A.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series in the matrix;
  --   pivots   pivots computed by the LU factorization.

  -- ON RETURN :
  --   Blead    leading exponents of the permuted matrix;
  --   Bcffs    coefficients of the matrix A, permuted with pivots.

  procedure Copy
              ( nrows,ncols,deg : in integer32;
                Alead : in Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : out Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Copies the matrix A into B.

  -- ON ENTRY :
  --   nrows    number of rows of the matrices;
  --   ncols    number of columns of the matrices;
  --   deg      degree of the series in the matrices;
  --   Alead    leading exponents of the series of the matrix A;
  --   Acffs    coefficients of the series in the matrix A;
  --   Blead    matrix of the same ranges as Alead;
  --   Bcffs    allocated space for the same dimensions as A.

  -- ON RETURN :
  --   Blead    copy of Alead;
  --   Bcffs    copy of Acffs.

  procedure Test_LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Tests the LU factorization with pivoting.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series in the matrix.

  procedure Test ( nrows,ncols,deg,low,upp : in integer32;
                   lower,upper : in boolean );

  -- DESCRIPTION :
  --   Generates a random linear system and runs a test.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   low      lower bound for the leading exponents;
  --   upp      upper bound for the leading exponents;
  --   lower    true if the test matrix is lower triangular;
  --   upper    true if the test matrix is upper triangular.

  procedure Specific_Test;

  -- DESCRIPTION :
  --   Runs a specific test.

  function Seed_Prompt return integer32;

  -- DESCRIPTION :
  --   Prompts for a fixed seed or a random seed.

  procedure Determinant_Test;

  -- DESCRIPTION :
  --   Tests the determinant of a random 2-by-2 matrix.

  procedure Random_Test;

  -- DESCRIPTION :
  --   Tests the LU factorization for a general random matrix.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the parameters of the tests and the runs tests.

end Test_Double_Lseries_Matrices;
