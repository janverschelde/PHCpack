with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Vectors;
with Double_Double_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;

package Test_DoblDobl_Linear_Solvers is

-- DESCRIPTION :
--   Tests solving linear systems of matrices of double doubles.
--   The second part of the code does a performance test.

  procedure Run_Double_Double_Linear_Solver
              ( A : in Double_Double_Matrices.Matrix;
                b : in Double_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Run_DoblDobl_Complex_Linear_Solver
              ( A : in DoblDobl_Complex_Matrices.Matrix;
                b : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test_Double_Double_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

  procedure Test_DoblDobl_Complex_Linear_Solver;

  -- DESCRIPTION :
  --   Generates a random linear system and solves it.

-- PERFORMANCE TESTERS :

  procedure Test_Performance_Standard_Float_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   This routine is primarily to compare how much more double double
  --   arithmetic costs, compared to standard hardware arithmetic.
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Test_Performance_Standard_Complex_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   This routine is primarily to compare how much more double double
  --   arithmetic costs, compared to standard hardware arithmetic.
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Test_Performance_Double_Double_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Test_Performance_DoblDobl_Complex_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean );

  -- DESCRIPTION :
  --   Generates a random matrix and random right hand side vector
  --   of dimension d, solves the system and reports the residual.
  --   This is done k times.
  --   If output, then the residual is computed and its max norm shown.

  procedure Ask_Performance_Settings
              ( n,d : out integer32;
                output,outer,solve : out boolean );

  -- DESCRIPTION :
  --   Prompts for the parameters of the performance tests:
  --   1) n : number of times to solve a linear system;
  --   2) d : dimension of the linear system to be solved;
  --   3) output : if the residual should be shown;
  --   4) outer : if a new random system should be generated each time;
  --   5) solve : solve or only LU.

  procedure Standard_Float_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems and times the solver.

  procedure Standard_Complex_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems and times the solver.

  procedure Double_Double_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems of double doubles
  --   and times the solver.

  procedure DoblDobl_Complex_Performance_Test;

  -- DESCRIPTION :
  --   Generates a number of linear systems
  --   and times the solver.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_DoblDobl_Linear_Solvers;
