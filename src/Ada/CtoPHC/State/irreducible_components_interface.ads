with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Irreducible_Components_Interface is

-- DESCRIPTION :
--   Provides the numerical irreducible decomposition, wrapping the
--   blackbox solver to run cascade homotopies, to filter and factor.

  function Standard_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in standard
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function Standard_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in standard
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function DoblDobl_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in double
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function DoblDobl_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in double
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function QuadDobl_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in quad
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function QuadDobl_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in quad
  --   double precision, eventually followed by a filter and a factor.

  -- ON ENTRY :
  --   a       in a[0] is the top dimension of the solution set,
  --           in a[1] is the number of tasks, 0 for no multitasking,
  --           in a[2] is 0 or 1, if witness supersets need filtering,
  --           in a[3] is 0 or 1, if witness sets need factoring;
  --   b       in b[0] is the verbose flag, 0 for false, 1 for true,
  --           if verbose, then output is written to screen;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       b[0] has the size of the string representation of
  --           the irreducible factors, if factoring was requested.

  function Standard_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   standard double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Standard_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   standard double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function DoblDobl_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   double double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function DoblDobl_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   double double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function QuadDobl_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   quad double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function QuadDobl_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   quad double precision and copies it into the containers.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the witness set;
  --   vrblvl  is the verbose level.

  function Standard_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates the witness solutions in standard double precision.
  --   The verbose level is given in vrblvl.

  function DoblDobl_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates the witness solutions in double double precision.
  --   The verbose level is given in vrblvl.

  function QuadDobl_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates the witness solutions in quad double precision.
  --   The verbose level is given in vrblvl.

  function Irreducible_Factor_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Given in b the number of characters in the string representation
  --   of the irreducible factors, returns this string representation.
  --   Assigns to a[0], the number of characters in the computed
  --   string representation, for later verification.
  --   The verbose level is given in vrblvl.

end Irreducible_Components_Interface;
