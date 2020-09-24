with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Solutions;
with TripDobl_Speelpenning_Convolutions;

package TripDobl_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius in triple double precision
--   for artificial or natural-parameter homotopies.

  procedure TripDobl_Newton_Fabry
              ( cfs : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in TripDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs.

  procedure TripDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out TripDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure TripDobl_Artificial_Setup;

  -- DESCRIPTION :
  --   Promps for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined.

  procedure TripDobl_Natural_Setup;

  -- DESCRIPTION :
  --   Promps for a natural-parameter homotopy, with start solutions.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the type of homotopy.

end TripDobl_Fabry_on_Homotopy;
