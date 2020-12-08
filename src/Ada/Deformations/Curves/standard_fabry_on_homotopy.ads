with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;

package Standard_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius in double precision
--   for artificial-parameter or natural-parameter homotopies.

  procedure Newton_Fabry
              ( cvh : Standard_Speelpenning_Convolutions.Link_to_System;
                cfs : in Standard_Coefficient_Convolutions.Link_to_System;
                sol : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs,
  --   on coefficient convolution systems (if no multitasking)
  --   with a staggered, inlined, and indexed Newton's method.
  --   With multitasking, the original convolution systems are used.

  procedure Run ( nbequ,idxpar,deg : in integer32;
                  sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With the homotopy defined starting at the solutions in sols,
  --   runs Newton's method on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

  procedure Run ( file : in file_type; nbequ,idxpar,deg : in integer32;
                  sols : in out Standard_Complex_Solutions.Solution_List );

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

end Standard_Fabry_on_Homotopy;
