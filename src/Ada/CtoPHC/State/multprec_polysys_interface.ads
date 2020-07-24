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

  function Multprec_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all polynomials with multiprecision coefficients.
  --   The verbose level is given in vrblvl.

end Multprec_PolySys_Interface;
