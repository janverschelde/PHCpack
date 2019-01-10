with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Series_Path_Trackers is

-- DESCRIPTION :
--   A series path tracker runs Newton's method on power series to
--   compute Pade approximants to predict solutions in the path tracker.
--   The procedures in this package give access to such trackers
--   in double, double double, and quad double precision.

  procedure Standard_Run
              ( nq,nvr,idxpar : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Run
              ( nq,nvr,idxpar : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Run
              ( nq,nvr,idxpar : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With a homotopy defined, runs the path tracker starting at the
  --   given solutions, in double, double double or quad double precision.

  -- ON ENTRY :
  --   nq       number of equations in the homotopy;
  --   nvr      number of variables in the homotopy;
  --   idxpar   index of the parameter in a natural parameter homotopy,
  --            in 1..nvr, or else 0 for an artificial parameter homotopy.

  function Prompt_for_Artificial return boolean;

  -- DESCRIPTION :
  --   Asks the user whether the homotopy is an artificial parameter
  --   homotopy and return true if so, otherwise false is returned.

  procedure Standard_Main;
  procedure DoblDobl_Main;
  procedure QuadDobl_Main;

  -- DESCRIPTION :
  --   Prompts the user for a homotopy, runs the series path trackers
  --   in standard double, double double, or quad double precision.

end Series_Path_Trackers;
