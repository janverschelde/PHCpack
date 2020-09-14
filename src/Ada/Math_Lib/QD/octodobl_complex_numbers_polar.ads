with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;

package OctoDobl_Complex_Numbers_Polar is

-- DESCRIPTION :
--   Provides a polar view on the octo double complex numbers.

  function Radius ( c : Complex_Number ) return octo_double;

  -- DESCRIPTION :
  --   Returns the radius of the complex number.

  function Angle ( c : Complex_Number ) return octo_double;

  -- DESCRIPTION :
  --   Returns the angle of the complex number.

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the ith root of the equation x^n - c = 0.

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number;
  function Polar_Exponentiation
             ( x : Complex_Number; e : octo_double ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x.

end OctoDobl_Complex_Numbers_Polar;
