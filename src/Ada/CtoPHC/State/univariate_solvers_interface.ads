with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Univariate_Solvers_Interface is

-- DESCRIPTION :
--   The functions below interface to the univariate root finders
--   in double, double double, quad double, and arbitrary precision.

  function Standard_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Calls the univariate polynomial solver in standard precision.
  --   The polynomial must be available in the container.
  --   The maximum number of iterations and the accuracy requirement
  --   are extracted from the input arguments a and c respectively.
  --   The solutions are placed in the solution container.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   c       in c[1] is the accuracy requirement;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of iterations done.

  function DoblDobl_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Calls the univariate polynomial solver in double double precision.
  --   The polynomial must be available in the container.
  --   The maximum number of iterations and the accuracy requirement
  --   are extracted from the input arguments a and c respectively.
  --   The solutions are placed in the solution container.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   c       in c[1] is the accuracy requirement;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of iterations done.

  function QuadDobl_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Calls the univariate polynomial solver in quad double precision.
  --   The polynomial must be available in the container.
  --   The maximum number of iterations and the accuracy requirement
  --   are extracted from the input arguments a and c respectively.
  --   The solutions are placed in the solution container.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   c       in c[1] is the accuracy requirement;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of iterations done.

  function Multprec_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Calls the univariate polynomial solver in multiprecision.
  --   The polynomial must be available in the container.
  --   The maximum number of iterations and the accuracy requirement
  --   are extracted from the input arguments a and c respectively.
  --   The solutions are placed in the solution container.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places in the precision,
  --           in a[1] is the maximum number of iterations;
  --   c       in c[1] is the accuracy requirement;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of iterations done.

end Univariate_Solvers_Interface;
