with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Scaling_Interface is

-- DESCRIPTION :
--   The functions in this package define the interface to the
--   equation and coefficient scaling of systems and solutions.

  function Scale_Standard_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in standard double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, 2 to define the type of scaling:
  --           0 : only equation scaling,
  --           1 : variable scaling without variability reduction
  --           2 : variable scaling with variability reduction;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       if a[0] is 1 or 2, then c contains the
  --           scaling coefficients and the estimate for the
  --           condition number, stored as a standard complex vector,
  --           the size of this vector is 4*n + 2.

  function Scale_DoblDobl_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, 2 to define the type of scaling:
  --           0 : only equation scaling,
  --           1 : variable scaling without variability reduction
  --           2 : variable scaling with variability reduction;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       if a[0] is 1 or 2, then c contains the
  --           scaling coefficients and the estimate for the
  --           condition number, stored as a dobldobl complex vector,
  --           the size of this vector is 8*n + 4.

  function Scale_QuadDobl_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, 2 to define the type of scaling:
  --           0 : only equation scaling,
  --           1 : variable scaling without variability reduction
  --           2 : variable scaling with variability reduction;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       if a[0] is 1 or 2, then c contains the
  --           scaling coefficients and the estimate for the
  --           condition number, stored as a quadobl complex vector,
  --           the size of this vector is 16*n + 8.

  function Scale_Multprec_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in arbitrary multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, 2 to define the type of scaling:
  --           0 : only equation scaling,
  --           1 : variable scaling without variability reduction
  --           2 : variable scaling with variability reduction;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       if a[0] is 1 or 2, then c contains the
  --           scaling coefficients and the estimate for the
  --           condition number, stored as a quadobl complex vector,
  --           the size of this vector is 16*n + 8.

  function Scale_Standard_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the standard solutions container.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the solution vector,
  --           but then as the dimension of all its doubles;
  --   b       defines the basis of the scaling;
  --   c       doubles for use in the coefficient scaling,
  --           note that a[0]/2 is the dimension of the complex vector
  --           and a[0]/2 is the dimension of the solution vector.

  function Scale_DoblDobl_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the dobldobl solutions container.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the solution vector,
  --           but then as the dimension of all its doubles;
  --   b       defines the basis of the scaling;
  --   c       doubles for use in the coefficient scaling,
  --           note that a[0]/4 is the dimension of the complex vector
  --           and a[0]/4 is the dimension of the solution vector.

  function Scale_QuadDobl_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the quaddobl solutions container.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the solution vector,
  --           but then as the dimension of all its doubles;
  --   b       defines the basis of the scaling;
  --   c       doubles for use in the coefficient scaling,
  --           note that a[0]/8 is the dimension of the complex vector
  --           and a[0]/8 is the dimension of the solution vector.

  function Scale_Multprec_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the multprec solutions container.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the solution vector,
  --           but then as the dimension of all its doubles;
  --   b       defines the basis of the scaling;
  --   c       doubles for use in the coefficient scaling,
  --           note that a[0]/8 is the dimension of the complex vector
  --           and a[0]/8 is the dimension of the solution vector.

end Scaling_Interface;
