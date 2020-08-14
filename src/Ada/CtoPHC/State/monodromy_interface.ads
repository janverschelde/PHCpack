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

end Monodromy_Interface;
