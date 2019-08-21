with text_io;                            use text_io;
with Ada.Calendar;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Homotopy_Continuation_Parameters;

package Drivers_to_Series_Trackers is

-- DESCRIPTION :
--   The procedures in this package help launching the path trackers
--   which apply power series methods as predictors.

  procedure Standard_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 );

  -- DESCRIPTION :
  --   Resets the gamma with a new Standard_Homotopy.Create.

  procedure DoblDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 );

  -- DESCRIPTION :
  --   Resets the gamma with a new DoblDobl_Homotopy.Create.

  procedure QuadDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 );

  -- DESCRIPTION :
  --   Resets the gamma with a new QuadDobl_Homotopy.Create.

  procedure Set_Output
              ( file : in out file_type;
                monitor,verbose,tofile : out boolean );

  -- DESCRIPTION :
  --   Prompts the user if verbose or not, and if so, whether the output
  --   should be written to a file or not, which sets tofile to true.
  --   If tofile, then file is created, ready for output.
  --   If the user wants to monitor the progress of the trackers,
  --   then the flag monitor is true on return.

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

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
  --   mhom     0 for affine coordinates, 1 for 1-homogenization,
  --            and m for m-homogenization;
  --   idz      index representation of the partition of m-homogenization,
  --            in case mhom > 1;
  --   verbose  if extra output during the tracking is needed;
  --   vrblvl   the verbose level.

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

  procedure Write_Conclusion 
              ( file : in file_type; start_moment : in Ada.Calendar.Time );

  -- DESCRIPTION :
  --   Given the moment of the start of the computations, writes the time
  --   stamp for the start and current moment and the seed number to file.

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the root refiners on the solutions in the list sols,
  --   where the nq is the number of equations.
  --   Output is written to the file.
  --   The verbose level is given in vrblvl.

  procedure Refine_Roots
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Refine_Roots
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Refine_Roots
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the root refiners on the solutions in the list sols,
  --   in double, double double, or quad double precision,
  --   converting the solutions to affine coordinates if mhom > 0.
  --   The mixed residuals are computed, evaluating the solutions
  --   at the homotopy polynomials with positive coefficients in abh.

  -- ON ENTRY :
  --   file     for output of the refined solutions;
  --   abh      the homotopy polynomials with absolute coefficients
  --            to compute the mixed residuals;
  --   mhom     0 if affine, 1 if 1-homogeneous, m if m-homogeneous;
  --   idz      the index representation of the partition, for mhom > 1;
  --   sols     solutions which will be refived;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sols     the refined solutions.

end Drivers_to_Series_Trackers;
