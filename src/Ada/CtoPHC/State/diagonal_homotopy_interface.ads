with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Diagonal_Homotopy_Interface is

-- DESCRIPTION :
--   The functions below give access to diagonal homotopies
--   to intersect solution sets.

  function Diagonal_Homotopy_Standard_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a polynomial in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Standard_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a Laurent polynomial
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a Laurent polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a polynomial in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a Laurent polynomial
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a Laurent polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a polynomial in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a witness set for a Laurent polynomial
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string;
  --   b       the string representation of a Laurent polynomial;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Standard_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for polynomial systems
  --   stored as target and start systems in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for polynomial systems
  --   stored as target and start systems in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for polynomial systems
  --   stored as target and start systems in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Standard_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for Laurent systems
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for Laurent systems
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a diagonal homotopy for Laurent systems
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first set;
  --   b       in b[0] is the dimension of the second set;
  --   vrblvl  is the verbose level.

end Diagonal_Homotopy_Interface;
