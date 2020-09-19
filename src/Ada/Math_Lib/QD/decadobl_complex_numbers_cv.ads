with Standard_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package DecaDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex deca doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_DecaDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_DecaDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a deca double complex number.

  function DecaDobl_Complex_to_Standard
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function DecaDobl_Complex_to_Multprec
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a deca double complex number
  --   into a standard or multiprecision complex number.

end DecaDobl_Complex_Numbers_cv;
