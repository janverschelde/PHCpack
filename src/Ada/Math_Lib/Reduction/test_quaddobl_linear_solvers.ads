with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Vectors;
with Quad_Double_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package Test_QuadDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests the LU decomposition on quad double matrices.

  procedure Run_Quad_Double_Linear_Solver
              ( A : in Quad_Double_Matrices.Matrix;
                b : in Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_QuadDobl_Complex_Linear_Solver
              ( A : in QuadDobl_Complex_Matrices.Matrix;
                b : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Quad_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_QuadDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_Performance_Quad_Double_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Test_Performance_QuadDobl_Complex_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Quad_Double_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems
  --   and times the solver.

  procedure QuadDobl_Complex_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems
  --   and times the solver.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_QuadDobl_Linear_Solvers;
