with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Drivers_to_Cascade_Filtering is

-- DESCRIPTION :
--   This package contains procedures to set up the sequence of
--   homotopies to compute generic points on all components.

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pocotime : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pocotime : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pocotime : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pocotime : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pocotime : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pocotime : out duration );

  -- DESCRIPTION :
  --   Performs a continuation to remove the slice from the embedded system.
  --   On entry, sols contains the start solutions, on return, the
  --   computed solutions are in the list sols,
  --   in standard double, double double, or quad double precision.
  --   The number of tasks equals nt.  For no tasking set nt = 0.

  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   resfile   file for final results:
  --               1) system at each level, for k down to 0; and
  --               2) solutions with zero slack variables;
  --   nt        number of tasks, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float );

  -- DESCRIPTION :
  --   This witness generate writes the witness supersets to files,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Standard_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );
  procedure DoblDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );
  procedure QuadDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Interactive driver to call the Witness_Generate procedure,
  --   in standard double, double double, or quad double precision.
  --   Names for the input and output file are passed at the command line.

  -- ON ENTRY :
  --   nt       number of tasks, if zero, then no tasking;
  --   inpname  name of the input file;
  --   outname  name of the output file.

  procedure Driver_to_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Prompts the user for the level of the working precision and
  --   then calls the Standard, DoblDobl, or QuadDobl_Witness_Generate.
  --   Names for the input and output file are passed at the command line.

  -- ON ENTRY :
  --   nt       number of tasks, if zero, then no tasking;
  --   inpname  name of the input file;
  --   outname  name of the output file.

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Laur_Systems.Laur_Sys );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Prompts the user to enter the top dimension, embeds the system,
  --   and then runs a cascade of homotopies,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     file opened for output;
  --   name     name of the output file;
  --   nt       number of tasks for multitasking,
  --            if zero, then no multitasking will be used;
  --   p        a polynomial system.

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Does the embedding of the top dimension and runs the cascade,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nt       the number of tasks, if 0 then no multitasking,
  --            otherwise nt tasks will be used to track the paths;
  --   inpname  name of the input file;
  --   outname  name of the output file.

  procedure Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Prompts the user for the level of working precision,
  --   does the embedding of the top dimension and runs the cascade.

  -- ON ENTRY :
  --   nt       the number of tasks, if 0 then no multitasking,
  --            otherwise nt tasks will be used to track the paths;
  --   inpname  name of the input file;
  --   outname  name of the output file.

end Drivers_to_Cascade_Filtering;
