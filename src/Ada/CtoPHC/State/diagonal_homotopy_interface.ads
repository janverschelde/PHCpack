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

  function Diagonal_Homotopy_Symbols_Doubler
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Doubles the number of symbols in the symbol table to write the
  --   solved target system to start the cascade of diagonal homotopies
  --   in extrinsic coordinates.
  --   On a succesful return, the symbol table contains the suffixed
  --   symbols to write the target system properly.

  -- ON ENTRY :
  --   a       in a[0] is the ambient dimension,
  --           equal to the original number of variables,
  --           in a[1] is the top dimension of the set,
  --           in a[2] is the number of characters in the string b;
  --   b       the sequence of symbols for the variables in the first set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Standard_Start_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions to start the cascade to intersect
  --   two witness sets in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first witness set;
  --   b       in b[0] is the dimension of the second witness set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Start_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions to start the cascade to intersect
  --   two witness sets in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first witness set;
  --   b       in b[0] is the dimension of the second witness set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Start_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions to start the cascade to intersect
  --   two witness sets in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the dimension of the first witness set;
  --   b       in b[0] is the dimension of the second witness set;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Prompt_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts for the first or second witness set to be read from file.

  -- ON ENTRY :
  --   a       in a[0] is the index of the set (one or two);
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the dimension of the witness set,
  --           in b[1] is the degree of the witness set.

  function Diagonal_Homotopy_Reset_Input
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Resets the input file for a witness set.

  -- ON ENTRY :
  --   a       in a[0] is the index of the witness set;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the degree of the witness set,
  --           in b[1] is the dimension of the witness set.

  function Diagonal_Homotopy_Cascade_Dimension
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to intersect two
  --   witness sets in extrinsic coordinates.

  -- ON ENTRY :
  --   a       in a[0] is the ambient dimension of the first set,
  --           in a[1] is the ambient dimension of the second set;
  --   b       in b[0] is the dimension of the first set,
  --           in b[1] is the dimension of the second set;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the dimension of the extrinsic cascade.

  function Diagonal_Homotopy_Standard_Hyperset
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   For one of the polynomials stored in double precision,
  --   makes a witness set for the hypersurface defined by that polynomial.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial in double precision,
  --           in a[1] is the number of characters in the string b;
  --   b       contains the name of the output file;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_Standard_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Elimates the extrinsic diagonal from system and solutions
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables,
  --           in a[1] is the number of slack variables to be added;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_DoblDobl_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Elimates the extrinsic diagonal from system and solutions
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables,
  --           in a[1] is the number of slack variables to be added;
  --   vrblvl  is the verbose level.

  function Diagonal_Homotopy_QuadDobl_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Elimates the extrinsic diagonal from system and solutions
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of slack variables,
  --           in a[1] is the number of slack variables to be added;
  --   vrblvl  is the verbose level.

end Diagonal_Homotopy_Interface;
