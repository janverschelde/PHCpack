with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Step_Trackers_Interface is
 
  function Step_Trackers_Standard_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the homotopy to track paths in double precision.

  -- ON ENTRY :
  --   a       if a[0] = 1, then a fixed gamma constant is used,
  --           otherwise the gamma is provided in c;
  --   c       in c[0] is the real part of the gamma constant,
  --           in c[1] is the imaginary part of the gamma constant;
  --   vrblvl  is the verbose level.

  function Step_Trackers_DoblDobl_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the homotopy to track paths in double double precision.

  -- ON ENTRY :
  --   a       if a[0] = 1, then a fixed gamma constant is used,
  --           otherwise the gamma is provided in c;
  --   c       in c[0] is the real part of the gamma constant,
  --           in c[1] is the imaginary part of the gamma constant;
  --   vrblvl  is the verbose level.

  function Step_Trackers_QuadDobl_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the homotopy to track paths in quad double precision.

  -- ON ENTRY :
  --   a       if a[0] = 1, then a fixed gamma constant is used,
  --           otherwise the gamma is provided in c;
  --   c       in c[0] is the real part of the gamma constant,
  --           in c[1] is the imaginary part of the gamma constant;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Set_Standard_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the start solution for the step-by-step tracker
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Set_DoblDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the start solution for the step-by-step tracker
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in double double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Set_QuadDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the start solution for the step-by-step tracker
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in quad double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Next_Standard_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next solution in double precision,
  --   replaces the current solution in double precision
  --   with the next solution.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Next_DoblDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next solution in double double precision,
  --   replaces the current solution in double double precision
  --   with the next solution.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in double double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Next_QuadDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next solution in quad double precision,
  --   replaces the current solution in quad double precision
  --   with the next solution.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in quad double precision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Multprec_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the homotopy to track paths in multiprecision.

  -- ON ENTRY :
  --   a       if a[0] = 1, then a fixed gamma constant is used;
  --   b       in b[0] is the number of decimal places;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Set_Multprec_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the start solution for the step-by-step tracker
  --   in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in multiprecision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Next_Multprec_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next solution in multiprecision,
  --   replaces the current solution in multiprecision
  --   with the next solution.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution stored
  --           in multiprecision;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Varbprec_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the homotopy to track paths in multiprecision.

  -- ON ENTRY :
  --   a       if a[0] = 1, then a fixed gamma constant is used,
  --           in a[1] is the total number of characters in the string b,
  --           in a[2] is the number of characters in the first string
  --           that contains the target system;
  --   b       the string representation of two systems, respectively
  --           the target and start system in the homotopy;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Set_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the start solution for the step-by-step tracker
  --   in variable multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in b,
  --           in a[1] is the number of variables in the solution,
  --           given in the string b;
  --   b       the string representation of a solution;
  --   vrblvl  is the verbose level.

  function Step_Trackers_Next_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next solution in variable multiprecision,
  --   replaces the current solution in variable multiprecision
  --   with the next solution.

  -- ON ENTRY :
  --   a       in a[0] is the wanted number of accurate decimal places,
  --           in a[1] is the maximum precision to be used,
  --           in a[2] is the maximum number of corrector steps,
  --           in a[3[ is 1 or 0, respectively corresponding to
  --           whether intermediate output is wanted or not.
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the size of the string representation of the next solution,
  --           to retrieve the string representation of the next solution,
  --           call Step_Trackers_Get_Varbprec_Solution.

  function Step_Trackers_Get_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the current variable precision solution as a string,
  --   as computed by Step_Trackers_Next_Varbprec_Solution.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in b;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the current solution;

  function Step_Trackers_Standard_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data allocated for the step-by-step trackers
  --   in double precision.

  function Step_Trackers_DoblDobl_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data allocated for the step-by-step trackers
  --   in double double precision.

  function Step_Trackers_QuadDobl_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data allocated for the step-by-step trackers
  --   in quad double precision.

  function Step_Trackers_Multprec_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data allocated for the step-by-step trackers
  --   in multiprecision.

  function Step_Trackers_Varbprec_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data allocated for the step-by-step trackers
  --   in variable multiprecision.

end Step_Trackers_Interface;
