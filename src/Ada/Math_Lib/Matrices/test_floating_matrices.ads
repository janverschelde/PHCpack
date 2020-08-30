package Test_Floating_Matrices is

-- DESCRIPTION :
--   Tests the matrix packages of standard and multi-precision floats.

  procedure Test_Standard_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of floating-point numbers,
  --   in standard double precision.  Writes the matrix to screen.

  procedure Test_Standard_VecMat_io;

  -- DESCRIPTION :
  --   Prompts for a number of floating-point matrices.
  --   Writes the vector of matrices to screen.

  procedure Test_Multprec_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then for the
  --   entries of a matrix of multiprecision floating-point numbers.
  --   Writes the matrix to screen.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Floating_Matrices;
