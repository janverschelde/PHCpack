with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;

package Test_Integer_Inverse is

-- DESCRIPTION :
--   Test procedures to compute the inverse of an integer matrix.

  procedure Test_Standard_Inverse
              ( A : in Standard_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Computes the inverse of the square matrix A.

  -- REQUIRED : A'range(1) = A'range(2).

  procedure Test_Standard_Inverse
              ( A : in Standard_Integer64_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Computes the inverse of the square matrix A.

  -- REQUIRED : A'range(1) = A'range(2).

  procedure Test_Multprec_Inverse
              ( A : in Multprec_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Computes the inverse of the square matrix A.

  -- REQUIRED : A'range(1) = A'range(2).

  function Standard_Unimodular_Smith_Transformation
              ( A : Standard_Integer_Matrices.Matrix ) 
              return Standard_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the Smith normal form of A with standard 32-bit arithmetic
  --   and returns the product of the unimodular transformations used to
  --   bring A in diagonal form.

  function Standard_Unimodular_Smith_Transformation
              ( A : Standard_Integer64_Matrices.Matrix ) 
              return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the Smith normal form of A with standard 64-bit arithmetic
  --   and returns the product of the unimodular transformations used to
  --   bring A in diagonal form.

  function Multprec_Unimodular_Smith_Transformation
              ( A : Multprec_Integer_Matrices.Matrix ) 
              return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the Smith normal form of A with multiprecision arithmetic
  --   and returns the product of the unimodular transformations used to 
  --   bring A in diagonal form.

  procedure Inverse_of_Standard32_Random_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for lower and uper bounds on the entries of a random
  --   n-dimensional matrix of 32-bit integers and runs a number of tests.

  procedure Inverse_of_Standard64_Random_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for lower and uper bounds on the entries of a random
  --   n-dimensional matrix of 64-bit integers and runs a number of tests.

  procedure Inverse_of_Standard_Random_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for 32-bit or 64-bit arithmetic
  --   and then runs a test on a random n-dimensional matrix.

  procedure Inverse_of_Multprec_Random_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for lower and uper bounds on the entries of a random
  --   n-dimensional matrix of multiprecision integers and then
  --   runs a number of tests.

  procedure Inverse_of_Given_Standard32_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for an n-by-n integer matrix and attempts
  --   to compute its inverse, using standard 32-bit arithmetic.

  procedure Inverse_of_Given_Standard64_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for an n-by-n integer matrix and attempts
  --   to compute its inverse, using 64-bit arithmetic.

  procedure Inverse_of_Given_Standard_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for an n-by-n integer matrix and attempts
  --   to compute its inverse.

  procedure Inverse_of_Given_Multprec_Matrix ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts for an n-by-n integer matrix and attempts
  --   to compute its inverse.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Integer_Inverse;
