with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package DoblDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex double doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_DoblDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_DoblDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a double double complex number.

  function DoblDobl_Complex_to_Standard
             ( c : DoblDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function DoblDobl_Complex_to_Multprec
             ( c : DoblDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a double double complex number
  --   into a standard or multiprecision complex number.

end DoblDobl_Complex_Numbers_cv;
