with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Multprec_PolySys_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   polynomials with multiprecision coefficients.
--   The integer returned by all functions should be zero if the job
--   was successful, otherwise the nonzero job code is returned.

  function Multprec_PolySys_Read
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with multiprecision
  --   coefficients and initializes the container with the system.
  --   The verbose level is given in vrblvl.

  function Multprec_PolySys_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a system with quad double precision coefficients
  --   from file into the systems container.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --           in a[1] is the working precision to parse the coefficients;
  --   b       the characters in the file name;
  --   vrblvl  the verbose level.

  function Multprec_PolySys_Write
             ( vrblvl : in integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the polynomial system with multiprecision coefficients
  --   to the defined output file or to screen.
  --   The verbose level is given in vrblvl.
 
  function Multprec_PolySys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the polynomial system
  --   with multiprecision coefficients.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of polynomials with multiprecision coefficients.

  function Multprec_PolySys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the number of polynomials with multiprecision
  --   coefficients and initializes the symbol table as well.

  -- ON ENTRY :
  --   a       the dimension of the polynomial system;
  --   vrblvl  the verbose level.

  function Multprec_PolySys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of terms of a polynomial
  --   with multiprecision coefficients.

  -- ON ENTRY :
  --   a       a[1] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       a[0] is the number of terms of a polynomial 
  --           stored in multiprecision with index in a[1].

  function Multprec_PolySys_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of a polynomial stored
  --   in multiprecision.

  -- ON ENTRY :
  --   a       a[0] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       b[0] is the degree of a polynomial 
  --           stored in multiprecision with index in a[0].

  function Multprec_PolySys_String_Save
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Given the string representation of a polynomial,
  --   stores the polynomial in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string,
  --           in a[1] is the number of variables in the polynomial,
  --           in a[2] is the index of the polynomial in the system;
  --   b       the string representation of the polynomial;
  --   vrblvl  is the verbose level.

  function Multprec_PolySys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size limit of the string representation of
  --   a polynomial in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the size limit of the string representation of the
  --           polynomial in multprecision with index a[0].

  function Multprec_PolySys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of the polynomial
  --   stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of characters in b;
  --   b       the string representation of the polynomial.

  function Multprec_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all polynomials with multiprecision coefficients.
  --   The verbose level is given in vrblvl.

  function Multprec_PolySys_Prompt_for_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a target system in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places in the precision;
  --   vrblvl  is the verbose level.

  function Multprec_PolySys_Write_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the target system in multiprecision.
  --   The verbose level is in vrblvl.

  function Multprec_PolySys_Prompt_for_Start
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a start system in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of decimal places in the precision;
  --   vrblvl  is the verbose level.

  function Multprec_PolySys_Write_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the start system in multiprecision.
  --   The verbose level is in vrblvl.

end Multprec_PolySys_Interface;
