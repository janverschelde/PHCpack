with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;
with DoblDobl_Speelpenning_Convolutions;

package DoblDobl_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius on double double precision
--   for artificial-parameter or natural-parameter homotopies.

  procedure Newton_Fabry
              ( nbt : in natural32;
                cfs : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs.
  --   In this interactive version, the user is prompted each time 
  --   for the parameters and all output is written to screen.
  --   The number of tasks is in nbt.

  procedure Run ( nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.
  --   All output is written to screen in this interactive run.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure Run ( file : in file_type;
                  nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  maxit : in integer32; tol : in double_float;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  verbose : in boolean );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.
  --   All output is written to file and the run is not interactive.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbt      the number of tasks;
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on the update vector;
  --   sols     start solutions in the homotopy.

  procedure Artificial_Setup
              ( nbtasks : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined.
  --   The number of tasks is in nbtasks
  --   and the verbose level in vrblvl.

  procedure Natural_Setup
              ( nbtasks : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for a natural-parameter homotopy, with start solutions.
  --   The number of tasks is in nbtasks
  --   and the verbose level in vrblvl.

  procedure Main ( nbtasks : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the type of homotopy.
  --   The number of tasks is in nbtasks
  --   and the verbose level in vrblvl.

end DoblDobl_Fabry_on_Homotopy;
