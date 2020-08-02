with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Reduction_Interface is

-- DESCRIPTION :
--   The functions below interface to the reduction methods to
--   reduce the total degree of the polynomial system.

  function Reduction_Standard_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system stored in double precision,
  --   and applies linear reduction in double precision.
  --   The system is replaced by the reduced system.

  -- ON ENTRY :
  --   a       if a[1] = 1, then the coefficient matrix will diagonalized,
  --           otherwise, the coefficient matrix will be made triangular;
  --   vrblvl  is the verbose level.

  function Reduction_DoblDobl_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system stored in double double precision,
  --   and applies linear reduction in double double precision.
  --   The system is replaced by the reduced system.

  -- ON ENTRY :
  --   a       if a[1] = 1, then the coefficient matrix will diagonalized,
  --           otherwise, the coefficient matrix will be made triangular;
  --   vrblvl  is the verbose level.

  function Reduction_QuadDobl_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system stored in double double precision,
  --   and applies linear reduction in double double precision.
  --   The system is replaced by the reduced system.

  -- ON ENTRY :
  --   a       if a[1] = 1, then the coefficient matrix will diagonalized,
  --           otherwise, the coefficient matrix will be made triangular;
  --   vrblvl  is the verbose level.

  function Reduction_Standard_Nonlinear
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --  Applies nonlinear reduction to the system in double precision.
  --  The system is replaced by the reduced system.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of equal degree replacements,
  --           in a[1] is the maximum number of computed S-polynomials, and
  --           in a[2]  sthe maximum number of computed R-polynomials;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of equal degree replacements,
  --           in b[1] is the number of computed S-polynomials, and
  --           in b[2] is the number of computed R-polynomials.

end Reduction_Interface;
