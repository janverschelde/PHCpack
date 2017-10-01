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

package Cascade_Homotopies is

-- DESCRIPTION :
--   A cascade homotopy defines a sequence of homotopies to compute
--   to compute generic points on all components of the solution set.
--   Versions of the cascade homotopies are provided
--   1) for ordinary and for Laurent polynomial systems, and
--   2) in double, double double, and quad double precision.
--   The six versions are combined with various levels of output:
--   1) the output is written to one single file, or
--   2) each superwitness set is written to a separate file.

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

end Cascade_Homotopies;
