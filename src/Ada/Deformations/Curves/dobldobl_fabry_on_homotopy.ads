with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;
with DoblDobl_Speelpenning_Convolutions;

package DoblDobl_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius on double double precision
--   for artificial or natural-parameter homotopies.

  procedure DoblDobl_Newton_Fabry
              ( cfs : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs.

  procedure DoblDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure DoblDobl_Artificial_Setup;

  -- DESCRIPTION :
  --   Promps for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined.

  procedure DoblDobl_Natural_Setup;

  -- DESCRIPTION :
  --   Promps for a natural-parameter homotopy, with start solutions.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the type of homotopy.

end DoblDobl_Fabry_on_Homotopy;
