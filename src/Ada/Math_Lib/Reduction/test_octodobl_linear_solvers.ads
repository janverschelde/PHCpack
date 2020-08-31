with Octo_Double_Vectors;
with Octo_Double_Matrices;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;

package Test_OctoDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of octo doubles.

  procedure Run_Octo_Double_Linear_Solver
              ( A : in Octo_Double_Matrices.Matrix;
                b : in Octo_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_OctoDobl_Complex_Linear_Solver
              ( A : in OctoDobl_Complex_Matrices.Matrix;
                b : in OctoDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Octo_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_OctoDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_OctoDobl_Linear_Solvers;
