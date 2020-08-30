package Test_Integer_Matrices is

-- DESCRIPTION :
--   Tests matrices with 32-bit, 64-bit, and multiprecision integers,
--   and vectors of integer matrices as well.

  procedure Test_Standard_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of 32-bit integer numbers.
  --   Writes the matrix to screen.

  procedure Test_Standard64_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of 64-bit integer numbers.
  --   Writes the matrix to screen.

  procedure Test_Standard_VecMat_io;

  -- DESCRIPTION :
  --   Prompts for a number of 32-bit integer matrices.
  --   Writes the vector of matrices to screen.

  procedure Test_Standard64_VecMat_io;

  -- DESCRIPTION :
  --   Prompts for a number of 64-bit integer matrices.
  --   Writes the vector of matrices to screen.

  procedure Test_Multprec_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of multiprecision integer numbers.
  --   Writes the matrix to screen.

  procedure Test_Multprec_Matrix_Vector_Product;

  -- DESCRIPTION :
  --   Prompts for a matrix and a vector and computes the
  --   matrix times vectors with multiprecision arithmetic.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Integer_Matrices;
