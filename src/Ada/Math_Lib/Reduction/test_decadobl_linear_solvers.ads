with Deca_Double_Vectors;
with Deca_Double_Matrices;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;

package Test_DecaDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of deca doubles.

  procedure Run_Deca_Double_Linear_Solver
              ( A : in Deca_Double_Matrices.Matrix;
                b : in Deca_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_DecaDobl_Complex_Linear_Solver
              ( A : in DecaDobl_Complex_Matrices.Matrix;
                b : in DecaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Deca_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_DecaDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_DecaDobl_Linear_Solvers;
