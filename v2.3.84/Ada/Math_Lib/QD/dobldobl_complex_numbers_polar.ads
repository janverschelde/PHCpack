with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;

package DoblDobl_Complex_Numbers_Polar is

-- DESCRIPTION :
--   Offers a polar view on the double double complex numbers.

  function Radius ( c : Complex_Number ) return double_double;

  -- DESCRIPTION :
  --   Returns the radius of the complex number.

  function Angle ( c : Complex_Number ) return double_double;

  -- DESCRIPTION :
  --   Returns the angle of the complex number.

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the ith root of the equation x^n - c = 0.

end DoblDobl_Complex_Numbers_Polar;
