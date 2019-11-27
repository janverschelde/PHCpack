with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Solutions;
with Path_Parameters;                    use Path_Parameters;

package Algorithmic_DiffEval_Trackers is

-- DESCRIPTION :
--   Defines the interface to the path trackers which use the reverse mode
--   of algorithmic differentiation.  There are two types of procedures:
--   (1) Wrappers to the C interface function to the Path library;
--   (2) Interactive drivers to test the wrappers.
--   Exported methods are Newton's method, a tracker for one path,
--   and the tracking of many paths. 
--   The level of precision is limited to standard double precision.

-- WRAPPERS :

  procedure Standard_ADE_Newton
              ( verbose : in integer32; pars : in Parameters );

  -- DESCRIPTION :
  --   Wraps the call to the C interface to run Newton's method in
  --   standard double, double double, or quad double precision
  --   on the target system stored in the systems container and the
  --   first solution stored in the solutions container.
  --   If verbose is positive, then output is written to screen.
  --   The values for the path parameters are in pars.

  -- REQUIRED :
  --   The containers have been initialized with a target system
  --   and at least one start solution.

  procedure Standard_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters );

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers, tracks one path
  --   in standard double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy;
  --   pars     values for the path parameters.

  -- REQUIRED :
  --   The containers have been initialized with a target system,
  --   a start system, and at least one start solution.

  procedure Standard_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters );

  -- DESCRIPTION :
  --   Wraps the C interface function to track many paths.
  --   With the homotopy as defined by the containers, tracks paths
  --   in standard double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy;
  --   pars     values for the path parameters.

  -- REQUIRED :
  --   The containers have been initialized with a target system,
  --   a start system, and a list of start solutions.

-- INTERACTIVE DRIVERS :

  procedure Standard_Newton;

  -- DESCRIPTION :
  --   This interactive procedure prompts the user for a system,
  --   and a corresponding list of solutions.
  --   Then the wrapper to the C interface Newton method is executed,
  --   in standard double precision.

  procedure Standard_Track_one_Path;

  -- DESCRIPTION :
  --   This interactive procedure prompts the user for a target system,
  --   start system, and a corresponding start solution.
  --   Then the wrapper to the C interface path tracker is executed,
  --   in standard double precision.
  --   Even if the user would give a whole list of start solutions,
  --   only the path starting at the first solution is tracked.

  procedure Standard_Refine_Roots
              ( file : in file_type;
                sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the root refiners and the computation of the condition tables
  --   on the list of solutions for the target system stored in the
  --   package PHCpack_Operations.

  procedure Standard_Track_many_Paths;

  -- DESCRIPTION :
  --   This interactive procedure prompts the user for a target system,
  --   start system, and a corresponding list of start solutions.
  --   Then the wrapper to the C interface path tracker is executed,
  --   in standard double precision.
  --   As many paths are tracked as there are start solutions.

end Algorithmic_DiffEval_Trackers;
