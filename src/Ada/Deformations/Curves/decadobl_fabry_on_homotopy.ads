with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Solutions;
with DecaDobl_Speelpenning_Convolutions;

package DecaDobl_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius in deca double precision
--   for artificial-parameter or natural-parameter homotopies.

  procedure Newton_Fabry
              ( cfs : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DecaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs.
  --   In this interactive version, the user is prompted each time 
  --   for the parameters and all output is written to screen.

  procedure Run ( nbequ,idxpar,deg : in integer32;
                  sols : in out DecaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.
  --   All output is written to screen in this interactive run.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure Run ( file : in file_type; nbequ,idxpar,deg : in integer32;
                  sols : in out DecaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.
  --   All output is written to file and the run is not interactive.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure Artificial_Setup;

  -- DESCRIPTION :
  --   Prompts for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined.

  procedure Natural_Setup;

  -- DESCRIPTION :
  --   Prompts for a natural-parameter homotopy, with start solutions.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the type of homotopy.

end DecaDobl_Fabry_on_Homotopy;
