with text_io;                            use text_io;
with Ada.Calendar;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Main_Homotopy_Continuation is

-- DESCRIPTION :
--   Defines the homotopy, with an artificial or a natural parameter,
--   calls the path trackers, and then the root refiners.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out Standard_Complex_Solutions.Solution_List;
                solsfile : in boolean; nbequ : in natural32 := 0 );

  -- DESCRIPTION :
  --   Calls the standard root refiners or sharpeners in case nbequ > 0.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out DoblDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean );

  -- DESCRIPTION :
  --   Calls the double double precision root refiners.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out QuadDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean );

  -- DESCRIPTION :
  --   Calls the quad double precision root refiners.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out Standard_Complex_Solutions.Solution_List;
               -- solsfile : in boolean;
                nbequ : in natural32 := 0 );

  -- DESCRIPTION :
  --   Calls the standard root refiners or sharpeners in case nbequ > 0.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out DoblDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean );

  -- DESCRIPTION :
  --   Calls the double double precision root refiners.

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out QuadDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean );

  -- DESCRIPTION :
  --   Calls the quad double precision root refiners.

  procedure Write_Solutions
             ( solsft : in file_type;
               stsols : in Standard_Complex_Solutions.Solution_List;
               ddsols : in DoblDobl_Complex_Solutions.Solution_List;
               qdsols : in QuadDobl_Complex_Solutions.Solution_List;
               mpsols : in Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the nonempty solution list to file solsft.

  procedure Write_Conclusion
              ( file : in file_type;
                start_moment,ended_moment : in Ada.Calendar.Time );

  -- DESCRIPTION :
  --   Writes to file the time stamps for the wall clock time elapsed
  --   between the start and end moment of the computations.
  --   Writes to file also the seed and version number.

  procedure Secant_Polynomial_Homotopy
              ( outfilename : in string;
                nbequ,nbvar,prclvl : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Creates the output file and reads start system and start solutions
  --   for an artificial parameter homotopy for polynomial systems.

  procedure Secant_Laurent_Homotopy
              ( outfilename : in string;
                nbequ,nbvar,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Creates the output file and reads start system and start solutions
  --   for an artificial parameter homotopy for Laurent systems.

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines a setup for a parameter homotopy in standard double precision.

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines a setup for a parameter homotopy in double double precision.

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines a setup for a parameter homotopy in quad double precision.

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Sets the parameters and runs the sweep homotopy
  --   in standard double precision.

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Sets the parameters and runs the sweep homotopy
  --   in double double precision.

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Sets the parameters and runs the sweep homotopy
  --   in quad double precision.

  procedure Multitasking_Secant_Homotopy
               ( outfilename : in string;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 ls : in Link_to_Array_of_Strings;
                 nt,nbequ,nbvar,prclvl : in natural32;
                 vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the multitasking path trackers on the system p.

  -- ON ENTRY :
  --   outfilename is the name of the output file;
  --   p         the system parsed to standard double precision;
  --   ls        the string representation of the system p;
  --   nt        the number of tasks;
  --   nbequ     number of equations;
  --   nbvar     number of variables;
  --   prclvl    the precision level;
  --   verbose   the verbose level.

  procedure Multitasking_Secant_Homotopy
               ( outfilename : in string;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 ls : in Link_to_Array_of_Strings;
                 nt,nbequ,nbvar,prclvl : in natural32;
                 vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the multitasking path trackers on the Laurent system p.

  -- ON ENTRY :
  --   outfilename is the name of the output file;
  --   p         the system parsed to standard double precision;
  --   ls        the string representation of the system p;
  --   nt        the number of tasks;
  --   nbequ     number of equations;
  --   nbvar     number of variables;
  --   prclvl    the precision level;
  --   verbose   the verbose level.

  procedure Parameter_or_Sweep_Homotopy
              ( inft : in out file_type; prclvl : in natural32;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the type of homotopy: parameter or sweep,
  --   and asks for the level of precision, but only if prclvl = 1.

  procedure Polynomial_Tracker
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Main procedure to call the path trackers
  --   on a regular polynomial system.

  procedure Standard_Laurent_Tracker
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Main procedure to call the path trackers
  --   on a Laurent polynomial system.

  procedure Main_Dispatch
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Decides whether the given system p is a genuine Laurent system,
  --   that is: if there are negative exponents, and then calls either
  --   the path trackers for Laurent systems (not as well developed),
  --   or the path trackers for regular polynomial systems.
  --   The general path trackers may run in higher precision and use
  --   the original formulation of the system as given in the strings ls.
  --   The verbose level is in vrb.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   prclvl : in natural32; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads in the system as a pointer to an array of strings
  --   and converts to a Laurent polynomial system with standard
  --   complex coefficients for a first dispatch.
  --   Defines what phc -p does.

end Main_Homotopy_Continuation;
