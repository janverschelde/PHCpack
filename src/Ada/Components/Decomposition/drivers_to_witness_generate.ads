with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Drivers_to_Witness_Generate is

-- DESCRIPTION :
--   The original code for a numerical irreducible decomposition,
--   with the witness generate and classify method, is here.

  procedure Timing_Summary ( file : in file_type;
                             roco,hoco,poco,total : in duration );

  -- DESCRIPTION :
  --   Displays the timing summary for the black-box solver.

  procedure Black_Box_Solver
              ( file : in file_type;
                sys : in Standard_Complex_Poly_Systems.Poly_Sys;
                deg : in boolean;
                sols : out Standard_Complex_Solutions.Solution_List;
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

  procedure Write_Generate_Summary
              ( file : in file_type;
                npaths,nsols0,nsols1,ndiv : in Standard_Natural_Vectors.Vector;
                timings : in Array_of_Duration );

  -- DESCRIPTION :
  --   Writes the summary for the cascade of homotopies.

  procedure Write ( file : in file_type; fp : in List; i : in integer32 );

  -- DESCRIPTION :
  --   Writes the i-th component of every element in the list.

  procedure Write_Banner
               ( file : in file_type; m : in natural32; sep : character );

  -- DESCRIPTION :
  --   Auxiliary to write the banner to file in the classify summary,
  --   as a sequence of m*6 characters sep.

  procedure Write_Classify_Summary
              ( file : in file_type; fp : in List;
                timings : in Array_of_Duration );

  -- DESCRIPTION :
  --   Writes the summary of the classification of witness points.

  procedure Driver_for_Cascade_Filter
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32 );

  -- DESCRIPTION :
  --   Performs a cascade of homotopies to find generic points on
  --   the solution set of the system p.  k is the top dimension.
  --   All results are written to the file.

end Drivers_to_Witness_Generate;
