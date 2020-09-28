with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Test_Standard_Solver is

-- DESCRIPTION :
--   Tests the black box solver in double precision.

  procedure Solve_to_File
              ( nt : in natural32; mvonly : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the name of the output file and solves p.

  -- ON ENTRY :
  --   nt       the number of tasks, 0 for no multitasking;
  --   mvonly   focus on mixed volumes and polyhedral homotopies;
  --   p        a square polynomial system;
  --   vrb      the verbose level.

  procedure Solve_without_File
              ( nt : in natural32; mvonly : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Solves p without output to file.

  -- ON ENTRY :
  --   nt       the number of tasks, 0 for no multitasking;
  --   mvonly   focus on mixed volumes and polyhedral homotopies;
  --   p        a square polynomial system;
  --   vrb      the verbose level.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial systems and asks for the setup:
  --   (1) the verbose level;
  --   (2) the number of tasks;
  --   (3) focus on mixed volumes and polyhedral homotopies;
  --   (4) output to file or not.

end Test_Standard_Solver;
