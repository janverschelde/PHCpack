with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Deflation_Interface is

-- DESCRIPTION :
--   The functions in this package interface to the deflation methods
--   to recondition isolated singular solutions.

  function Deflation_Standard_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies deflation in double precision to the polynomial system
  --   and to the solutions stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   b       in b[0] is the maximum number of deflations;
  --   c       in c[0] is the tolerance on the forward error,
  --           in c[1] is the tolerance on the backward error,
  --           in c[2] is the tolerance on the rank;
  --   vrblvl  is the verbose level.

  function Deflation_DoblDobl_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies deflation in double double precision to the polynomial
  --   system and to the solutions stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   b       in b[0] is the maximum number of deflations;
  --   c       in c[0] is the tolerance on the forward error,
  --           in c[1] is the tolerance on the backward error,
  --           in c[2] is the tolerance on the rank;
  --   vrblvl  is the verbose level.

  function Deflation_QuadDobl_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Applies deflation in double double precision to the polynomial
  --   system and to the solutions stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum number of iterations;
  --   b       in b[0] is the maximum number of deflations;
  --   c       in c[0] is the tolerance on the forward error,
  --           in c[1] is the tolerance on the backward error,
  --           in c[2] is the tolerance on the rank;
  --   vrblvl  is the verbose level.

  function Deflation_Standard_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the multiplicity structure for the polynomial system
  --   and the solutions stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum differentiation order;
  --   b       in b[0] is the value for the verbose flag;
  --   c       in c[0] is the tolerance for the numerical rank;
  --   vrblvl  is the verbose level.

  -- ON RETRUN :
  --   a       in a[0] is the computed multiplicity;
  --   b       contains the multiplicity structure.

  function Deflation_DoblDobl_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the multiplicity structure for the polynomial system
  --   and the solutions stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum differentiation order;
  --   b       in b[0] is the value for the verbose flag;
  --   c       in c[0] is the tolerance for the numerical rank;
  --   vrblvl  is the verbose level.

  -- ON RETRUN :
  --   a       in a[0] is the computed multiplicity;
  --   b       contains the multiplicity structure.

  function Deflation_QuadDobl_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the multiplicity structure for the polynomial system
  --   and the solutions stored in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the maximum differentiation order;
  --   b       in b[0] is the value for the verbose flag;
  --   c       in c[0] is the tolerance for the numerical rank;
  --   vrblvl  is the verbose level.

  -- ON RETRUN :
  --   a       in a[0] is the computed multiplicity;
  --   b       contains the multiplicity structure.

end Deflation_Interface;
