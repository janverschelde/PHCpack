with Penta_Double_Vectors;
with Penta_Double_Matrices;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Matrices;

package Test_PentDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of penta doubles.

  procedure Run_Penta_Double_Linear_Solver
              ( A : in Penta_Double_Matrices.Matrix;
                b : in Penta_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_PentDobl_Complex_Linear_Solver
              ( A : in PentDobl_Complex_Matrices.Matrix;
                b : in PentDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Penta_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_PentDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_PentDobl_Linear_Solvers;
