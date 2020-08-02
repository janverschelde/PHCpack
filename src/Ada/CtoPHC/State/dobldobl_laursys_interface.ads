with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package DoblDobl_LaurSys_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   Laurent polynomials in double double precision.
--   The integer returned by all functions should be zero if the job
--   was successful, otherwise the nonzero job code is returned.

  function DoblDobl_LaurSys_Read
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system in double double
  --   precision and initializes the container with the system provided
  --   by the user.  The verbose level is given in vrblvl.

  function DoblDobl_LaurSys_Write
             ( vrblvl : in integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the polynomial system stored in double precision
  --   to the defined output file or to screen.
  --   The verbose level is given in vrblvl.
 
  function DoblDobl_LaurSys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the Laurent polynomial system stored
  --   in double double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of Laurent polynomials
  --           stored in double double precision.

  function DoblDobl_LaurSys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the number of Laurent polynomials in double double
  --   precision and initializes the symbol table as well.

  -- ON ENTRY :
  --   a       the dimension of the Laurent polynomial system;
  --   vrblvl  the verbose level.

  function DoblDobl_LaurSys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of terms of a Laurent polynomial stored
  --   in double double precision.

  -- ON ENTRY :
  --   a       a[1] is the index of a Laurent polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       a[0] is the number of terms of a Laurent polynomial 
  --           stored in double precision with index in a[1].

  function DoblDobl_LaurSys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a term of a Laurent polynomial 
  --   stored in double double precision.

  -- ON ENTRY :
  --   a       a[1] contains the number of the Laurent polynomial,
  --           a[2] is the index to a term in that Laurent polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       contains the exponent vector of the term;
  --   c       real and imaginary part of the coefficient.

  function DoblDobl_LaurSys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds a term to a Laurent polynomial stored in double double precision.

  -- ON ENTRY :
  --   a       a[1] contains the index of the Laurent polynomial;
  --   b       the exponents in the term;
  --   c       real and imaginary part of the coefficient;
  --   vrblvl  the verbose level.

  function DoblDobl_LaurSys_String_Save
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Given the string representation of a Laurent polynomial,
  --   stores the Laurent polynomial in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string,
  --           in a[1] is the number of variables in the polynomial,
  --           in a[2] is the index of the polynomial in the system;
  --   b       the string representation of the Laurent polynomial;
  --   vrblvl  is the verbose level.

  function DoblDobl_LaurSys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size limit of the string representation of
  --   a Laurent polynomial in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the Laurent polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the size limit of the string representation of the Laurent
  --           polynomial in double double precision with index a[0].

  function DoblDobl_LaurSys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of a Laurent polynomial
  --   stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of characters in b;
  --   b       the string representation of a polynomial.

  function DoblDobl_LaurSys_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops a variable from the Laurent system
  --   stored in double double precision.
  --   The variable to be dropped is given by its index.

  -- ON ENTRY :
  --   a      in a[0] is the index of the variable;
  --   vrblvl is the verbose level.

  function DoblDobl_LaurSys_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops a variable from the Laurent system 
  --   stored in double double precision.
  --   The variable to be dropped is given by its name.

  -- ON ENTRY :
  --   a      in a[0] is the number of characters in b;
  --   b      contains the name of the variable to be dropped;
  --   vrblvl is the verbose level.

  function DoblDobl_LaurSys_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all Laurent polynomials stored in double double precision.
  --   The verbose level is given in vrblvl.

end DoblDobl_LaurSys_Interface;
