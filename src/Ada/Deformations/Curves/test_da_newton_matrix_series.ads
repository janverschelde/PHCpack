package Test_DA_Newton_Matrix_Series is

-- DESCRIPTION :
--   Tests the application of Newton's method to compute series solutions,
--   in deca double precision.

  procedure Test_LU_Newton;

  -- DESCRIPTION :
  --   Prompts for a system of series coefficients
  --   and for an initial series approximation.
  --   Then, runs the tests on the LU Newton method.

  procedure Test_QR_Newton;

  -- DESCRIPTION :
  --   Prompts for a system of series coefficients
  --   and for an initial series approximation.
  --   Then, runs the tests on the QR Newton method.

  procedure Test_SVD_Newton;

  -- DESCRIPTION :
  --   Prompts for a system of series coefficients
  --   and for an initial series approximation.
  --   Then, runs the tests on the SVD Newton method.

  procedure Test_Echelon_Newton;

  -- DESCRIPTION :
  --   Prompts for a system of series coefficients
  --   and for an initial series approximation.
  --   Then, runs the tests on the Echelon Newton method.

  procedure Main;

  -- DESCRIPTION :
  --   Displays the test menu and then runs the selected test.

end Test_DA_Newton_Matrix_Series;
