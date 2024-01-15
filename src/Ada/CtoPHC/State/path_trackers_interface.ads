with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Path_Trackers_Interface is

-- DESCRIPTION :
--   Provides functions to define homotopies and run path trackers.

  function Path_Trackers_Standard_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double precision with random gamma constant.
  --   The verbose level is in vrblvl.

  function Path_Trackers_Standard_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double precision with given gamma constant.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Standard_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the homotopy in double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_DoblDobl_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double double precision with random gamma constant.
  --   The verbose level is given in verbose level.

  function Path_Trackers_DoblDobl_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double double precision with given gamma constant.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the homotopy in double double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_QuadDobl_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in quad double precision with given gamma constant.
  --   The verbose level is given in vrblvl.

  function Path_Trackers_QuadDobl_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in quad double precision with given gamma constant.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the homotopy in quad double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_Multprec_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in multiprecision with random gamma constant.
  --   The verbose level is given in vrblvl.

  function Path_Trackers_Multprec_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in multiprecision with given gamma constant.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Multprec_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the homotopy in multiprecision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_Standard_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double precision with random gamma constant,
  --   to solve a Laurent system.  The verbose level is in vrblvl.

  function Path_Trackers_Standard_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double precision with given gamma constant,
  --   for Laurent systems.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Standard_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the Laurent homotopy in double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_DoblDobl_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double double precision with random gamma constant,
  --   to solve a Laurent system.  The verbose level is in vrblvl.

  function Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in double double precision with given gamma constant,
  --   for Laurent systems.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the Laurent homotopy in double double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_QuadDobl_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in quad double precision with random gamma constant,
  --   to solve a Laurent system.  The verbose level is in vrblvl.

  function Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a homotopy in quad double precision with given gamma constant,
  --   for Laurent systems.

  -- ON ENTRY :
  --   a       in a[0] is the power of t in the homotopy;
  --   c       in c[0] is the real part of the gamma,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data for the Laurent homotopy in quad double precision.
  --   The verbose level is in vrblvl.

  function Path_Trackers_Standard_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in double precision without intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved;
  --   b       multiplicity flag;
  --   c       coordinates of the solution at the end.

  function Path_Trackers_DoblDobl_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in double double precision without intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved;
  --   b       multiplicity flag;
  --   c       coordinates of the solution at the end.

  function Path_Trackers_QuadDobl_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in quad double precision without intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved;
  --   b       multiplicity flag;
  --   c       coordinates of the solution at the end.

  function Path_Trackers_Standard_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in double precision with intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved.

  function Path_Trackers_DoblDobl_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in double double precision with intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved;
  --   b       multiplicity flag;
  --   c       coordinates of the solution at the end.

  function Path_Trackers_QuadDobl_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path in quad double precision with intermediate output,
  --   starting at the given solution.

  -- ON ENTRY :
  --   b       multiplicity flag;
  --   c       coordinates of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of steps,
  --           in a[1] is the number of failures,
  --           in a[2] is the number of iterations,
  --           in a[3] is the number of linear systems solved;
  --   b       multiplicity flag;
  --   c       coordinates of the solution at the end.

  function Path_Trackers_Standard_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solution in double precision
  --   given with diagnostics to the defined output file.

  -- ON ENTRY :
  --   a       diagnostics of the path tracker;
  --   b       multiplicity flag;
  --   c       coordinates of the solution;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solution in double double precision
  --   given with diagnostics to the defined output file.

  -- ON ENTRY :
  --   a       diagnostics of the path tracker;
  --   b       multiplicity flag;
  --   c       coordinates of the solution;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solution in quad double precision
  --   given with diagnostics to the defined output file.

  -- ON ENTRY :
  --   a       diagnostics of the path tracker;
  --   b       multiplicity flag;
  --   c       coordinates of the solution;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Standard_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the polynomial homotopy defined in double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the polynomial homotopy defined in double double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the polynomial homotopy defined in quad double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Multprec_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the polynomial homotopy defined in multiprecision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Standard_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the Laurent homotopy defined in double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the Laurent homotopy defined in double double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   With the Laurent homotopy defined in quad double precision,
  --   runs the path trackers, with the number of tasks on input.

  -- ON ENTRY :
  --   a       in a[0] is the number of tasks;
  --   vrblvl  is the verbose level.

  function Path_Trackers_Standard_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the crude path tracker in double precision
  --   on the defined homotopy and stored start solutions.
  --   A crude path trackers does no postprocessing on the solutions.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  function Path_Trackers_DoblDobl_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the crude path tracker in double double precision
  --   on the defined homotopy and stored start solutions.
  --   A crude path trackers does no postprocessing on the solutions.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  function Path_Trackers_QuadDobl_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the crude path tracker in quad double precision
  --   on the defined homotopy and stored start solutions.
  --   A crude path trackers does no postprocessing on the solutions.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   vrblvl  is the verbose level.

end Path_Trackers_Interface;
