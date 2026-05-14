with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Double_Puiseux_Interface is

-- DESCRIPTION :
--   Defines functions to solve linear systems of real powered series
--   and to run Newton steps on real powered Laurent homotopies,
--   in double precision.
--   Required is that the Double/DCMPLX_VecVecs containers are well defined
--   and that there is a corresponding Laurent system in the container for
--   Laurent systems in double precision.

  procedure Run_Linear_Solver
              ( dim,nbr : in integer32; c : in C_dblarrs.Pointer;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the first nbr terms of the solution of a linear system
  --   of real powered series, of dimension dim, after the powers and
  --   the coefficients of the input system are extracted.
  --   The coefficients of the solution series are returned in c.

  function Linear_Solver
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes leading terms of the solution of a linear system
  --   of real powered series.

  -- ON ENTRY :
  --   a       number of terms of the solution desired is in a[0];
  --   c       has enough space for all terms of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coefficients of the solution series.

  function Newton_Steps
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Does a number of Newton steps to compute the solution
  --   of a real powered Laurent homotopy.

  -- ON ENTRY :
  --   a       number of steps is in a[0];
  --   c       has enough space for all terms of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coefficients of the solution series.

end Double_Puiseux_Interface;
