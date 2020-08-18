with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Monodromy_Interface is

-- DESCRIPTION :
--   Provides functions to apply monodromy.

  function Monodromy_Standard_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with system and solutions
  --   stored in double precision.
  
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with system and solutions
  --   stored in double double precision.
  
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with system and solutions
  --   stored in double double precision.
  
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the value of a coefficient in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   b       in b[0] is the coefficient in the equation;
  --   c       in c[0] is the real part of the coefficient,
  --           in c[1] is the imaginary part of the coefficient;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the value of a coefficient in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   b       in b[0] is the coefficient in the equation;
  --   c       in (c[0], c[1]) is the real part of the coefficient,
  --           in (c[2], c[3]) is the imaginary part of the coefficient;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the value of a coefficient in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   b       in b[0] is the coefficient in the equation;
  --   c       in c[0..3] is the real part of the coefficient,
  --           in c[4..7] is the imaginary part of the coefficient;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the gamma in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   c       in c[0] is the real part of gamma,
  --           in c[1] is the imaginary part of gamma;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the gamma in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   c       in c[0..1] is the real part of gamma,
  --           in c[2..3] is the imaginary part of gamma;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the gamma in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the linear equation;
  --   c       in c[0..3] is the real part of gamma,
  --           in c[4..7] is the imaginary part of gamma;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Sample
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new witness set in double precision,
  --   for the new coefficients of the linear equations.

  function Monodromy_DoblDobl_Sample
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new witness set in double double precision,
  --   for the new coefficients of the linear equations.

  function Monodromy_QuadDobl_Sample
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new witness set in quad double precision,
  --   for the new coefficients of the linear equations.

  function Monodromy_Standard_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps slices and solutions in double precision
  --   to turn back in the monodromy loop.

  function Monodromy_DoblDobl_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps slices and solutions in double double precision
  --   to turn back in the monodromy loop.

  function Monodromy_QuadDobl_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Swaps slices and solutions in quad double precision
  --   to turn back in the monodromy loop.

  function Monodromy_Standard_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded system in double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_DoblDobl_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded system in double double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_QuadDobl_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded system in quad double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_Standard_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the first solution list in double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_DoblDobl_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the first solution list in double double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_QuadDobl_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the first solution list in quad double precision
  --   from the sampler to the container.
  --   The verbose level is given by vrblvl.

  function Monodromy_Standard_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies solutions in double precision from the grid
  --   into the container.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution list;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies solutions in double double precision from the grid
  --   into the container.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution list;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies solutions in quad double precision from the grid
  --   into the container.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution list;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the permutations to process monodromy loops
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of monodromy loops;
  --   b       in b[0] is the degree of the solution set,
  --           in b[1] is the dimension of the solution set;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the permutations to process monodromy loops
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of monodromy loops;
  --   b       in b[0] is the degree of the solution set,
  --           in b[1] is the dimension of the solution set;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the permutations to process monodromy loops
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of monodromy loops;
  --   b       in b[0] is the degree of the solution set,
  --           in b[1] is the dimension of the solution set;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in double precision
  --   to the permutations for processing.
  --   The verbose level is given in vrblvl.

  function Monodromy_DoblDobl_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in double double precision
  --   to the permutations for processing.
  --   The verbose level is given in vrblvl.

  function Monodromy_QuadDobl_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in double double precision
  --   to the permutations for processing.
  --   The verbose level is given in vrblvl.

  function Monodromy_Standard_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new permutation by the last solution list
  --   in double precision copied to the permutations data.

  -- ON ENTRY:
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the computed permutation.

  function Monodromy_DoblDobl_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new permutation by the last solution list
  --   in double double precision copied to the permutations data.

  -- ON ENTRY:
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the computed permutation.

  function Monodromy_QuadDobl_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes a new permutation by the last solution list
  --   in quad double precision copied to the permutations data.

  -- ON ENTRY:
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the computed permutation.

  function Monodromy_Standard_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Updates the the decomposition with a new permutation,
  --   computed from a solution list in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the permutation;
  --   b       the permutation;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the previous number of groupings,
  --           in a[1] is the new number of groupings.

  function Monodromy_DoblDobl_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Updates the the decomposition with a new permutation,
  --   computed from a solution list in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the permutation;
  --   b       the permutation;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the previous number of groupings,
  --           in a[1] is the new number of groupings.

  function Monodromy_QuadDobl_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Updates the the decomposition with a new permutation,
  --   computed from a solution list in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the permutation;
  --   b       the permutation;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the previous number of groupings,
  --           in a[1] is the new number of groupings.

  function Monodromy_Standard_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the current decomposition computed with solution lists
  --   in double precision, to the defined output file or to screen.
  --   The verbose level is given in vrblvl.

  function Monodromy_DoblDobl_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the current decomposition computed with solution lists in
  --   double double precision, to the defined output file or to screen.
  --   The verbose level is given in vrblvl.

  function Monodromy_QuadDobl_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the current decomposition computed with solution lists in
  --   quad double precision, to the defined output file or to screen.
  --   The verbose level is given in vrblvl.

  function Monodromy_Standard_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

  function Monodromy_DoblDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

  function Monodromy_QuadDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

  function Monodromy_Standard_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the diagnostics of the grid in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the maximum error of the samples,
  --           in c[1] is the minimum distance between the samples.

  function Monodromy_DoblDobl_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the diagnostics of the grid in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the maximum error of the samples,
  --           in c[1] is the minimum distance between the samples.

  function Monodromy_QuadDobl_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the diagnostics of the grid in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the maximum error of the samples,
  --           in c[1] is the minimum distance between the samples.

  function Monodromy_Standard_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the difference between the trace and the actual sum
  --   for a factor in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of points in the factor;
  --   b       identifies the points in the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the difference between the trace and actual sum.

  function Monodromy_DoblDobl_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the difference between the trace and the actual sum
  --   for a factor in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of points in the factor;
  --   b       identifies the points in the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the difference between the trace and actual sum.

  function Monodromy_QuadDobl_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the difference between the trace and the actual sum
  --   for a factor in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of points in the factor;
  --   b       identifies the points in the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the difference between the trace and actual sum.

  function Monodromy_Standard_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the slices in the sampler in double precision
  --   with the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number to initialize the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the slices in the sampler in double double precision
  --   with the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number to initialize the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the slices in the sampler in quad double precision
  --   with the given number.

  -- ON ENTRY :
  --   a       in a[0] is the number to initialize the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Finds the index of a solution in a double precision slice.

  -- ON ENTRY :
  --   a       in a[0] is the label of the solution,
  --           in a[1] is the index of the slice;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       is the index of the solution in the slice.

  function Monodromy_DoblDobl_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Finds the index of a solution in a double double precision slice.

  -- ON ENTRY :
  --   a       in a[0] is the label of the solution,
  --           in a[1] is the index of the slice;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       is the index of the solution in the slice.

  function Monodromy_QuadDobl_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Finds the index of a solution in a quad double precision slice.

  -- ON ENTRY :
  --   a       in a[0] is the label of the solution,
  --           in a[1] is the index of the slice;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       is the index of the solution in the slice.

  function Monodromy_Standard_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the target slices in double precision
  --   to the slice stored with the given index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the target slices in double double precision
  --   to the slice stored with the given index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the target slices in quad double precision
  --   to the slice stored with the given index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one loop in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the start slice,
  --           in a[1] is the index of the target slice;
  --   b       in b[0] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the index of the target solution.

  function Monodromy_DoblDobl_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one loop in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the start slice,
  --           in a[1] is the index of the target slice;
  --   b       in b[0] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the index of the target solution.

  function Monodromy_QuadDobl_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one loop in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the start slice,
  --           in a[1] is the index of the target slice;
  --   b       in b[0] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the index of the target solution.

  function Monodromy_Standard_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of factors in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of factors.

  function Monodromy_DoblDobl_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of factors in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of factors.

  function Monodromy_QuadDobl_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of factors in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of factors.

  function Monodromy_Standard_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a factor in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       defines the points in the factor.

  function Monodromy_DoblDobl_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a factor in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       defines the points in the factor.

  function Monodromy_QuadDobl_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a factor in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the factor;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       defines the points in the factor.

  function Monodromy_Standard_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in double precision to silent.
  --   The verbose level is given in vrblvl.

  function Monodromy_DoblDobl_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in double double precision to silent.
  --   The verbose level is given in vrblvl.

  function Monodromy_QuadDobl_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in quad double precision to silent.
  --   The verbose level is given in vrblvl.

  function Monodromy_Standard_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in double precision to verbose.
  --   The verbose level is given in vrblvl.

  function Monodromy_DoblDobl_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in double double precision to verbose.
  --   The verbose level is given in vrblvl.

  function Monodromy_QuadDobl_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the state of monodromy in double double precision to verbose.
  --   The verbose level is given in vrblvl.

  function Monodromy_Standard_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a random complex numbers in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       contains real and imaginary parts of a double
  --           complex number.

  function Monodromy_DoblDobl_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a random complex numbers in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       contains real and imaginary parts of a double double
  --           complex number.

  function Monodromy_QuadDobl_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a random complex numbers in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       contains real and imaginary parts of a quad double
  --           complex number.

  function Monodromy_Standard_Add_Slice
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds new slices in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the total number of doubles,
  --           in a[1] is the dimension of the solution set,
  --           in a[2] is the ambient dimension;
  --   c       the coefficients of the slices;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Get_Slice
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the coefficients of the slices in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the total number of doubles,
  --           in a[1] is the dimension of the solution set,
  --           in a[2] is the ambient dimension;
  --   b       in b[0] is the index of the slice;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       the coefficients of the slices.

  function Monodromy_Standard_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with Laurent system and solutions
  --   stored in double precision.
 
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_DoblDobl_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with Laurent system and solutions
  --   stored in double double precision.
 
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_QuadDobl_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the sampler with Laurent system and solutions
  --   stored in quad double precision.
 
  -- ON ENTRY :
  --   a       in a[0] is the dimension of the set;
  --   vrblvl  is the verbose level.

  function Monodromy_Standard_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded Laurent system from the sampler
  --   to the container in double precision.
  --   The verbose level is given in vrblvl.

  function Monodromy_DoblDobl_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded Laurent system from the sampler
  --   to the container in double double precision.
  --   The verbose level is given in vrblvl.

  function Monodromy_QuadDobl_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the embedded Laurent system from the sampler
  --   to the container in quad double precision.
  --   The verbose level is given in vrblvl.

end Monodromy_Interface;
