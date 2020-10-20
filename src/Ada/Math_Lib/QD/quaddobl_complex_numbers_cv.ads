with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package QuadDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex double doubles,
--   complex triple doubles, and double and multiprecision complex vectors.

  function Standard_to_QuadDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_QuadDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a quad double complex number.

  function QuadDobl_Complex_to_Standard
             ( c : QuadDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function QuadDobl_Complex_to_DoblDobl
             ( c : QuadDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number;
  function QuadDobl_Complex_to_TripDobl
             ( c : QuadDobl_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number;
  function QuadDobl_Complex_to_Multprec
             ( c : QuadDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a quad double complex number into a double,
  --   double double, triple double, or multiprecision complex number.

end QuadDobl_Complex_Numbers_cv;
