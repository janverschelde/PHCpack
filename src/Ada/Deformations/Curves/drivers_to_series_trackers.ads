with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Homotopy_Continuation_Parameters;

package Drivers_to_Series_Trackers is

-- DESCRIPTION :
--   The procedures in this package help launching the path trackers
--   which apply power series methods as predictors.

  procedure Standard_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Resets the gamma with a new Standard_Homotopy.Create.

  procedure DoblDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Resets the gamma with a new DoblDobl_Homotopy.Create.

  procedure QuadDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Resets the gamma with a new QuadDobl_Homotopy.Create.

  procedure Set_Output
              ( file : in out file_type; verbose,tofile : out boolean );

  -- DESCRIPTION :
  --   Prompts the user if verbose or not, and if so, whether the output
  --   should be written to a file or not, which sets tofile to true.
  --   If tofile, then file is created, ready for output.

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in boolean := false );
  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false );
  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );
  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false );
  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Given the homotopy defined in Standard_Homotopy, or DoblDob_Homotopy,
  --   or QuadDobl_Homotopy, and the homotopy continuation parameters tuned,
  --   the series trackers start at the solutions in sols,
  --   in double, or double double, or quad double precision.

  -- ON ENTRY :
  --   file     optional file for writing statistics;
  --   nq       number of equations in the homotopy;
  --   sols     start solutions;
  --   pars     values for the homotopy continuation parameters,
  --            if omitted, then default values are used;
  --   verbose  if extra output during the tracking is needed.

  -- ON RETURN :
  --   sols     the solutions at the end of the tracking.

  procedure Write_Timer
              ( file : in file_type;
                numdeg,dendeg,precision : in natural32;
                timer : in Timing_Widget );

  -- DESCRIPTION :
  --   Writes the times with as well the degrees of numerator
  --   and denominator of the Pade approximants.
  --   The precision is 0, 1, or 2, respectively
  --   for double, double double, or quad double precision.

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

end Drivers_to_Series_Trackers;
