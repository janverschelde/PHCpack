with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Schubert_Interface is

-- DESCRIPTION :
--   Functions in this package define the interface to the 
--   numerical Schubert calculus to resolve intersection conditions
--   and to run Littlewood-Richardson homotopies.

  function Schubert_Intersection_Conditions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Resolves the intersection conditions for a Schubert problem.

  -- ON ENTRY :
  --   a       in a[0] is the ambient dimension,
  --           in a[1] is the dimension of the solution planes,
  --           in a[2] is the number of intersection conditions,
  --           in a[3] is the verbose flag, 0 for silent, 1 for output;
  --   b       contains the integers for the brackets,
  --           as many as a[1]*a[2];
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       the formal root count.

  function Standard_LR_Homotopies
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Littlewood-Richardson homotopies in double precision.

  function DoblDobl_LR_Homotopies
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Littlewood-Richardson homotopies in double double precision.

  function QuadDobl_LR_Homotopies
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs Littlewood-Richardson homotopies in quad double precision.

end Schubert_Interface; 
