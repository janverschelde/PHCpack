with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;

package TripDobl_Complex_Numbers_Polar is

-- DESCRIPTION :
--   Provides a polar view on the triple double complex numbers.

  function Radius ( c : Complex_Number ) return triple_double;

  -- DESCRIPTION :
  --   Returns the radius of the complex number.

  function Angle ( c : Complex_Number ) return triple_double;

  -- DESCRIPTION :
  --   Returns the angle of the complex number.

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the ith root of the equation x^n - c = 0.

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number;
  function Polar_Exponentiation
             ( x : Complex_Number; e : triple_double ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x.

end TripDobl_Complex_Numbers_Polar;
