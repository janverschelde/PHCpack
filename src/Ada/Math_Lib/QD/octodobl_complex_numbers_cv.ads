with Standard_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package OctoDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex octo doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_OctoDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_OctoDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a octo double complex number.

  function OctoDobl_Complex_to_Standard
             ( c : OctoDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function OctoDobl_Complex_to_Multprec
             ( c : OctoDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a octo double complex number
  --   into a standard or multiprecision complex number.

end OctoDobl_Complex_Numbers_cv;
