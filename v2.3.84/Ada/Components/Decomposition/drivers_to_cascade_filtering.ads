with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_to_Cascade_Filtering is

-- DESCRIPTION :
--   This package contains procedures to set up the sequence of
--   homotopies to compute generic points on all components.

  procedure Driver_to_Square_and_Embed;

  -- DESCRIPTION :
  --   A system is made square by
  --      adding extra hyperplanes, if it is underdetermined; or
  --      adding slack variables, if it is overdetermined.
  --   The embedding to capture k dimensional components is executed
  --   by adding k random hyperplanes and k slack variables.
  --   This interactive routine reads in a system and creates a new
  --   square and embedded polynomial system.

  procedure Driver_to_Remove_Embedding;

  -- DESCRIPTION :
  --   Removes embed and slack variables from the system read on input.
  --   This operation undoes the squaring and embedding.

  procedure Witness_Generate
               ( outfile,resfile : in file_type; ep : in Poly_Sys;
                 sols : in Solution_List; k : in natural32;
                 zerotol : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   resfile   file for final results:
  --               1) system at each level, for k down to 0; and
  --               2) solutions with zero slack variables;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   k         number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 ep : in Poly_Sys; sols : in Solution_List;
                 k : in natural32; zerotol : in double_float );

  -- DESCRIPTION :
  --   This witness generate writes the witness supersets to files.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   k         number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   zerotol   tolerance to decide whether a number is zero or not.

  procedure Driver_to_Witness_Generate;

  -- DESCRIPTION :
  --   Interactive driver to call the Witness_Generate procedure.

  procedure Black_Box_Solver
               ( file : in file_type; sys : in Poly_Sys;
                 deg : in boolean; sols : out Solution_List;
                 rc : out natural32; totaltime : out duration );
 
  -- DESCRIPTION :
  --   The black-box solver consists of the following stages:
  --     1) root-counting, when deg is true then degrees are used
  --     2) construction of start system; and
  --     3) path following to the actual target system.
  --   Several diagnostics and intermediate results are written to file.

  -- ON ENTRY :
  --   file      output file for diagnostics and intermediate results;
  --   sys       a square polynomial system;
  --   deg       flag to indicate when degrees or polytopes are used
  --             to count roots and construct the start system;
 
  -- ON RETURN :
  --   rc        root count equals the number of paths followed;
  --   sols      solutions found at the end of the paths;
  --   totaltime is the total cpu time used in computing the results.

  procedure Driver_for_Cascade_Filter
               ( file : in file_type; p : in Poly_Sys; k : in integer32 );
  
  -- DESCRIPTION :
  --   Performs a cascade of homotopies to find generic points on
  --   the solution set of the system p.  k is the top dimension.
  --   All results are written to the file.

  procedure Embed_and_Cascade;

  -- DESCRIPTION :
  --   Does the embedding of the top dimension and runs the cascade.

end Drivers_to_Cascade_Filtering;
