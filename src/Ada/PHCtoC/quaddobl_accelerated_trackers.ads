with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package QuadDobl_Accelerated_Trackers is

-- DESCRIPTION :
--   Wraps the C++ code to accelerate Newton's method and path tracking
--   in quad double precision.  The Ada procedures are in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure QuadDobl_GPU_Newton ( execmode,verbose : in integer32 );

  -- DESCRIPTION :
  --   Calls the accelerated Newton's method in quad double precision.
 
  -- REQUIRED :
  --   The containers for polynomial systems and for solutions
  --   in quad double precision contain a system and a solution.
  --   The computed solution is place in the solutions container.

  -- ON ENTRY :
  --   execmode   0 (both cpu and gpu), 1 (cpu only), or 2 (gpu only);
  --   verbose    if > 0, then additional output is written to screen.

  procedure QuadDobl_GPU_Track_One
                ( execmode,verbose : in integer32; gamma : in Complex_Number );

  -- DESCRIPTION :
  --   Calls the accelerated path tracker in quad double precision,
  --   for tracking one solution path.

  -- REQUIRED :
  --   A target and start system are defined in PHCpack_Operations
  --   and the solutions container contains one start solution.

  -- ON ENTRY :
  --   execmode   0 (both cpu and gpu), 1 (cpu only), or 2 (gpu only);
  --   verbose    if > 0, then additional output is written to screen;
  --   gamma      value of the random gamma constant.

  procedure QuadDobl_GPU_Track_Many
                ( execmode,verbose : in integer32; gamma : in Complex_Number );

  -- DESCRIPTION :
  --   Calls the accelerated path tracker in quad double precision,
  --   for tracking many solution paths.

  -- REQUIRED :
  --   A target and start system are defined in PHCpack_Operations
  --   and the start solutions are stored in the solutions container.

  -- ON ENTRY :
  --   execmode   0 (both cpu and gpu), 1 (cpu only), or 2 (gpu only);
  --   verbose    if > 0, then additional output is written to screen;
  --   gamma      value of the random gamma constant.

end QuadDobl_Accelerated_Trackers; 
