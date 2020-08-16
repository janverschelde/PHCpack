with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Newton_Interface is

-- DESCRIPTION :
--   The functions below interface to Newton's method on polynomial
--   and Laurent systems in double, double double, quad double, and
--   arbitrary precision arithmetic.

  function Newton_Standard_Polynomial_Verify
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Verifies the solutions of the polynomial system in double precision,
  --   with default settings of parameters, without deflation.
  --   If the output file is defined, then output is written to file.
  --   The verbose level is given on input.

  function Newton_Standard_Polynomial_Refine
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Newton's method to refine the given solution
  --   on the polynomial system stored in double precision.

  -- ON ENTRY :
  --   b       in b[0] is the multiplicity flag;
  --   c       stores the coordinates of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the new value of the multiplicity flag;
  --   c       updated coordinates of the solution.

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

  function Newton_Varbprec_Step
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does one Newton step on the system given as string,
  --   starting at the solution stored in multiprecision.
  --   The system may have negative exponents.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the polynomial system,
  --           in a[1] is the number of characters in the string b,
  --           in a[2] is the wanted number of accurate decimal places,
  --           in a[3] is the maximum number of iterations,
  --           in a[4] is the maximum number of decimal places
  --           in the working precision;
  --   b       contains the string representation of a polynomial system;
  --   vrblvl  is the verbose level.

  function Newton_Standard_SysPool_Refine
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Refines a system in the pool in double precision
  --   with the corresponding solutions in the pool in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the pool;
  --   vrblvl  is the verbose level.

end Newton_Interface;
