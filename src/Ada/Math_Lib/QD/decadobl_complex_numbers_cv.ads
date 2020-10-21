with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
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

  function DecaDobl_Complex_to_DoblDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function DecaDobl_Complex_to_TripDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number;
  function DecaDobl_Complex_to_QuadDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number;
  function DecaDobl_Complex_to_PentDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number;
  function DecaDobl_Complex_to_OctoDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the most significant parts of c,
  --   for double double, triple double, quad double, penta double,
  --   or octo double precision.

end DecaDobl_Complex_Numbers_cv;
