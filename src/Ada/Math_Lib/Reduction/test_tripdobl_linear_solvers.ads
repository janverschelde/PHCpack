with Triple_Double_Vectors;
with Triple_Double_Matrices;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;

package Test_TripDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of triple doubles.

  procedure Run_Triple_Double_Linear_Solver
              ( A : in Triple_Double_Matrices.Matrix;
                b : in Triple_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_TripDobl_Complex_Linear_Solver
              ( A : in TripDobl_Complex_Matrices.Matrix;
                b : in TripDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Triple_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_TripDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_TripDobl_Linear_Solvers;
