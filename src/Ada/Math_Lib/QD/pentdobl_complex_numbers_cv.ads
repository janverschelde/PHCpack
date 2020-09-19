with Standard_Complex_Numbers;
with PentDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package PentDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex penta doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_PentDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_PentDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a penta double complex number.

  function PentDobl_Complex_to_Standard
             ( c : PentDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function PentDobl_Complex_to_Multprec
             ( c : PentDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a penta double complex number
  --   into a standard or multiprecision complex number.

end PentDobl_Complex_Numbers_cv;
