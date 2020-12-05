with Ada.Calendar;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Laur_Systems;

package Mixed_Volume_Calculator is

-- DESCRIPTION :
--   Mixed volumes and mixed cell configurations define start systems
--   in polyhedral homotopies to solve random coefficient systems.
--   Several lifting strategies and mixed volume calculators are available.

  procedure Read_System
              ( filename : in string;
                lq : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   Attempts to open the file with the given name for reading of
  --   a polynomial system.
  
  function Prompt_for_Lifting return natural32;

  -- DESCRIPTION :
  --   Prompts the user to select a lifting strategy and returns this
  --   choice as a natural number between 0 and 5.

  procedure Call_MixedVol
               ( file : in file_type; nt : in natural32;
                 lq : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                 v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Asks for the precision and then calls then MixedVol.

  procedure Lift_Set_and_Run
              ( nt : in natural32; outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                lq : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                v : in integer32 := 0 );

  -- DESCRIPTION :
  --   Settles on a lifting, computes mixed volumes,
  --   and runs polyhedral homotopies if so desired.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -m.

end Mixed_Volume_Calculator;
