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
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, runs in standard double precision.

  procedure DoblDobl_Run
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, runs in double double precision.

  procedure QuadDobl_Run
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, runs in quad double precision.

  procedure Standard_Main;

  -- DESCRIPTION :
  --   Prompts the user for a homotopy with series coefficients,
  --   runs the series path trackers in standard double precision.

  procedure DoblDobl_Main;

  -- DESCRIPTION :
  --   Prompts the user for a homotopy with series coefficients,
  --   runs the series path trackers in double double precision.

  procedure QuadDobl_Main;

  -- DESCRIPTION :
  --   Prompts the user for a homotopy with series coefficients,
  --   runs the series path trackers in quad double precision.

end Series_Path_Trackers;
