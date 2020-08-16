with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Cascade_Homotopy_Interface is

-- DESCRIPTION :
--   Provides function to the cascade homotopies.

  function Cascade_Homotopy_Standard_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade polynomial homotopy in double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_DoblDobl_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade polynomial homotopy in double double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_QuadDobl_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade polynomial homotopy in quad double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_Standard_Laurent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade Laurent homotopy in double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_DoblDobl_Laurent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade Laurent homotopy in double double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_QuadDobl_Laurent
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines a cascade Laurent homotopy in quad double precision.
  --   The verbose level is given in vrblvl.

  function Cascade_Homotopy_Cut_Slack
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Removes the last slack variable in the embedded polynomial system
  --   stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the current number of slack variables,
  --           equal to the dimension of the witness set;
  --   vrblvl  is the verbose level.

end Cascade_Homotopy_Interface;
