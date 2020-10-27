with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;

package Test_Stable_Mixed_Volumes is

-- DESCRIPTION :
--   Interactive development of the computation of stable mixed volumes.

  procedure Check_Minimal_Degrees
              ( p : in Laur_Sys; is_Laur : out boolean );

  -- DESCRIPTION :
  --   Reports the minimal degrees of all polynomials in the system p.
  --   The is_Laur on return is true if p has negative exponents.

  procedure Check_Supports ( p : in Laur_Sys; is_done : out boolean );

  -- DESCRIPTION :
  --   Checks whether the origin belongs to each support.

  procedure Check_Zero_Types 
             ( n,r : in integer32; p : in Laur_Sys;
               mix : in Standard_Integer_Vectors.Link_to_Vector;
               mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Writes the zero types for all cells to screen.

  procedure Check_Zero_Types ( p : in Laur_Sys );

  -- DESCRIPTION :
  --   Reads a mixed subdivision from file and checks the zero types
  --   of all cells.

  procedure Call_Polyhedral_Continuation ( p : in Laur_Sys );

  -- DESCRIPTION :
  --   Prompts for file to write the random coefficient to
  --   and a regular mixed-cell configuration.
  --   Runs the polyhedral homotopies to compute all solutions
  --   of a random coefficient system for p.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Stable_Mixed_Volumes;
