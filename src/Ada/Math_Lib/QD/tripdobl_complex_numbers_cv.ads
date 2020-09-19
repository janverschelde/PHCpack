with Standard_Complex_Numbers;
with TripDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package TripDobl_Complex_Numbers_cv is

-- DESCRIPTION :
--   Converts between vectors of complex triple doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_TripDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number;
  function Multprec_to_TripDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision complex number
  --   into a triple double complex number.

  function TripDobl_Complex_to_Standard
             ( c : TripDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number;
  function TripDobl_Complex_to_Multprec
             ( c : TripDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts a triple double complex number
  --   into a standard or multiprecision complex number.

end TripDobl_Complex_Numbers_cv;
