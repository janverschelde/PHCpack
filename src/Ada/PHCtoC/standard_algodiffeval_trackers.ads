with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_AlgoDiffEval_Trackers is

-- DESCRIPTION :
--   Wraps the C++ code to apply algorithmic differentiation for evaluation
--   in Newton's method and path tracking in standard double precision.
--   The Ada procedures are in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure Standard_ADE_Newton ( verbose : in integer32 );

  -- DESCRIPTION :
  --   Calls Newton's method in standard double precision.
 
  -- REQUIRED :
  --   The containers for polynomial systems and for solutions
  --   in standard double precision contain a system and a solution.
  --   The computed solution is place in the solutions container.

  -- ON ENTRY :
  --   verbose    if > 0, then additional output is written to screen.

  procedure Standard_ADE_Track_One
                ( verbose : in integer32; gamma : in Complex_Number );

  -- DESCRIPTION :
  --   Calls the path tracker in standard double precision,
  --   for tracking one solution path.

  -- REQUIRED :
  --   A target and start system are defined in PHCpack_Operations
  --   and the solutions container contains one start solution.

  -- ON ENTRY :
  --   verbose    if > 0, then additional output is written to screen;
  --   gamma      value for the random gamma constant.

  procedure Standard_ADE_Track_Many
                ( verbose : in integer32; gamma : in Complex_Number );

  -- DESCRIPTION :
  --   Calls the path tracker in standard double precision,
  --   for tracking many solution paths.

  -- REQUIRED :
  --   A target and start system are defined in PHCpack_Operations
  --   and the start solutions are stored in the solutions container.

  -- ON ENTRY :
  --   verbose    if > 0, then additional output is written to screen;
  --   gamma      value for the random gamma constant.

end Standard_AlgoDiffEval_Trackers; 
