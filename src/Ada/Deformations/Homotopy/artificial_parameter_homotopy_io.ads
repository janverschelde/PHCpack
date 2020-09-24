with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Artificial_Parameter_Homotopy_io is

-- DESCRIPTION :
--   Provides procedures to prompt the user to provide file names
--   for the target, start system, and start solutions in an 
--   artificial-parameter homotopy, in double, double double,
--   triple double, and quad double precision.

  procedure get ( target : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out Standard_Complex_Solutions.Solution_List );
  procedure get ( target : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure get ( target : out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out TripDobl_Complex_Solutions.Solution_List );
  procedure get ( target : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  start : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                  sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a target system and for a start system with
  --   corresponding start solutions.

end Artificial_Parameter_Homotopy_io;
