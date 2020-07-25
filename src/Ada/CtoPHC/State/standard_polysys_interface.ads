with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Standard_PolySys_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   polynomials in double precision.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function Standard_PolySys_Read
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system in double precision and
  --   initializes the container with the system provided by the user.
  --   The verbose level is given in vrblvl.

  function Standard_PolySys_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a system with double precision coefficients
  --   from file into the systems container.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --   b       the characters in the file name;
  --   vrblvl  the verbose level.

  function Standard_PolySys_Write
             ( vrblvl : in integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the polynomial system stored in double precision
  --   to the defined output file or to screen.
  --   The verbose level is given in vrblvl.
 
  function Standard_PolySys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the polynomial system stored
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of polynomials stored in double precision.

  function Standard_PolySys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the number of polynomials in double precision
  --   and initializes the symbol table as well.

  -- ON ENTRY :
  --   a       the dimension of the polynomial system;
  --   vrblvl  the verbose level.

  function Standard_PolySys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of terms of a polynomial stored
  --   in double precision.

  -- ON ENTRY :
  --   a       a[1] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       a[0] is the number of terms of a polynomial 
  --           stored in double precision with index in a[1].

  function Standard_PolySys_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of a polynomial stored
  --   in double precision.

  -- ON ENTRY :
  --   a       a[0] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       b[0] is the degree of a polynomial 
  --           stored in double precision with index in a[0].

  function Standard_PolySys_Total_Degree
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the total degree of the system.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the product of the degrees of the polynomials
  --           with double precision coefficients.

  function Standard_PolySys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a term of a polynomial stored in double precision.

  -- ON ENTRY :
  --   a       a[1] contains the number of the polynomial,
  --           a[2] is the index to a term in that polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       contains the exponent vector of the term;
  --   c       real and imaginary part of the coefficient.

  function Standard_PolySys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds a term to a polynomial stored in double precision.

  -- ON ENTRY :
  --   a       a[1] contains the index of the polynomial;
  --   b       the exponents in the term;
  --   c       real and imaginary part of the coefficient;
  --   vrblvl  the verbose level.

  function Standard_PolySys_Clear_Symbols
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the symbol table.

  function Standard_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all polynomials stored in double precision.
  --   The verbose level is given in vrblvl.

end Standard_PolySys_Interface;
