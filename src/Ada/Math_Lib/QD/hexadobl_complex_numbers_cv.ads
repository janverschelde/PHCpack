with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package HexaDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between complex hexa doubles,
--   double and multiprecision complex numbers,
--   and other multiple double complex number types.

  function Standard_to_HexaDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return HexaDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_HexaDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return HexaDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a double or multiprecision complex number
  --   into a hexa double complex number.

  function HexaDobl_Complex_to_Standard
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_Multprec
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a hexa double complex number
  --   into a double or multiprecision complex number.

  function HexaDobl_Complex_to_DoblDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_TripDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_QuadDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_PentDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_OctoDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number;
  function HexaDobl_Complex_to_DecaDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the most significant parts of c,
  --   as a double double, triple double, quad double, penta double,
  --   octo double, or deca double complex number.

end HexaDobl_Complex_Numbers_cv;
