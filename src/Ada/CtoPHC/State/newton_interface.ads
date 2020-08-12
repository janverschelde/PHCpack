with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Newton_Interface is

-- DESCRIPTION :
--   The functions below interface to Newton's method on polynomial
--   and Laurent systems in double, double double, quad double, and
--   arbitrary precision arithmetic.

  function Newton_Standard_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in double precision.
  --   The verbose level is given on input.

  function Newton_DoblDobl_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in double double precision.
  --   The verbose level is given on input.

  function Newton_QuadDobl_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in quad double precision.
  --   The verbose level is given on input.

  function Newton_Multprec_Polynomial_Step
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places
  --           for the working precision;
  --   vrblvl  is the verbose level.

  function Newton_Standard_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the Laurent polynomial system and solutions
  --   stored in double precision.
  --   The verbose level is given on input.

  function Newton_DoblDobl_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in double double precision.
  --   The verbose level is given on input.

  function Newton_QuadDobl_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the Laurent polynomial system and solutions
  --   stored in quad double precision.
  --   The verbose level is given on input.

  function Newton_Multprec_Laurent_Step
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the polynomial system and solutions
  --   stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places
  --           in the working precision;
  --   vrblvl  is the verbose level.

end Newton_Interface;
