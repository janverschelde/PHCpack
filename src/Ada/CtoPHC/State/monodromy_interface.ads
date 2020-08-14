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

end Monodromy_Interface;
