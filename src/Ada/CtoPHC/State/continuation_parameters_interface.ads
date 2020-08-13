with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Continuation_Parameters_Interface is

-- DESCRIPTION :
--   Provides interface functions to the tuning of continuation parameters,
--   either interactively or via the setting of numerical values.

  function Continuation_Parameters_Ask_Values
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for the values of the parameters,
  --   writing to the defined output file, or to standard output.
  --   The verbose level is given by vrblvl.

  function Continuation_Parameters_Ask_Output_Level
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for the output level,
  --   writing to the defined output file, or to standard output.
  --   The verbose level is given by vrblvl.

  function Continuation_Parameters_Get_All
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves all values of the continuation parameters.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       is an array of 34 values of the continuation parameters,
  --           in the order as displayed by Continuation_Parameters_Show.

  function Continuation_Parameters_Set_All
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Given 34 values for the continuation parameters,
  --   sets all values of the continuation parameters.

  -- ON ENTRY :
  --   c       34 values for the continuation parameters,
  --           in the order as displayed by Continuation_Parameters_Show;
  --   vrblvl  is the verbose level.

  function Continuation_Parameters_Autotune
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the values of the parameter depending on the given level
  --   and the given number of decimal places.

  -- ON ENTRY :
  --   a       in a[0] is the level of expected difficulty of the paths,
  --           at a higher level, the path tracker will spend more work;
  --   b       in b[0] is the number of decimal places, which determines
  --           the corrector tolerances;
  --   vrblvl  is the verbose level.

  function Continuation_Parameters_Show
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Shows the values of all continuation parameters on standard output.
  --   The verbose level is given by the value of vrblvl.

  function Continuation_Parameters_Get_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the value of a continuation parameter given its index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the parameter, a number in 1..34;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the value for the parameter defined by the index.

  function Continuation_Parameters_Set_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets a value of continuation based on the input parameters.

  -- ON ENTRY :
  --   a       in a[0] is the index of the parameter, a number in 1..34;
  --   c       in c[0] is the value for the parameter defined by the index;
  --   vrblvl  is the verbose level.

end Continuation_Parameters_Interface;
