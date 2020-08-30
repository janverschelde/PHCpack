package Test_Complex_Matrices is

-- DESCRIPTION :
--   Tests matrices of double precision and multiprecision complex numbers.

  procedure Test_Standard_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of complex numbers,
  --   in standard double precision.
  --   Writes the matrix to screen.

  procedure Test_Standard_VecMat_io;

  -- DESCRIPTION :
  --   Prompts for a number of complex matrices in double precision.
  --   Writes the matrices to screen.

  procedure Test_Multprec_io;

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns and then
  --   for the entries of a matrix of multiprecision complex numbers.
  --   Writes the matrix to screen.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Complex_Matrices;
