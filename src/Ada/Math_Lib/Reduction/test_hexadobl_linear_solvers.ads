with Hexa_Double_Vectors;
with Hexa_Double_Matrices;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Matrices;

package Test_HexaDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of hexa doubles.

  procedure Run_Hexa_Double_Linear_Solver
              ( A : in Hexa_Double_Matrices.Matrix;
                b : in Hexa_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_HexaDobl_Complex_Linear_Solver
              ( A : in HexaDobl_Complex_Matrices.Matrix;
                b : in HexaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Hexa_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_HexaDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_HexaDobl_Linear_Solvers;
