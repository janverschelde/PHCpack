package Crude_Path_Trackers is

-- DESCRIPTION :
--   A crude path tracker does not refine or do any kind of postprocessing
--   of the solutions at the end of the paths.
--   Solutions are directly appended to the solutions container,
--   which allows to monitor the progress of the path tracking.

  procedure Standard_Track_Paths ( verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks all paths starting at the solutions
  --   stored as start solutions, in double precision.
  --   The solutions container will contain all solutions.
  --   Solutions are appended directly after the computation.
  --   There is no root postprocessing stage.

  -- ON ENTRY :
  --   verbose  if true, then each solution vector is written to screen,
  --            otherwise, the path tracker remains mute.

  -- REQUIRED :
  --   The data in PHCpack_Operations have been initialized with
  --   a target system, start system with start solutions,
  --   in standard double precision.

  procedure DoblDobl_Track_Paths ( verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks all paths starting at the solutions
  --   stored as start solutions, in double double precision.
  --   The solutions container will contain all solutions.
  --   Solutions are appended directly after the computation.
  --   There is no root postprocessing stage.

  -- ON ENTRY :
  --   verbose  if true, then each solution vector is written to screen,
  --            otherwise, the path tracker remains mute.

  -- REQUIRED :
  --   The data in PHCpack_Operations have been initialized with
  --   a target system, start system with start solutions,
  --   in double double precision.

  procedure QuadDobl_Track_Paths ( verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks all paths starting at the solutions
  --   stored as start solutions, in quad double precision.
  --   The solutions container will contain all solutions.
  --   Solutions are appended directly after the computation.
  --   There is no root postprocessing stage.

  -- ON ENTRY :
  --   verbose  if true, then each solution vector is written to screen,
  --            otherwise, the path tracker remains mute.

  -- REQUIRED :
  --   The data in PHCpack_Operations have been initialized with
  --   a target system, start system with start solutions,
  --   in quad double precision.

end Crude_Path_Trackers;
