with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Complex_Numbers_Polar is

-- DESCRIPTION :
--   Offers a polar view on the standard complex numbers.

  function Radius ( c : Complex_Number ) return double_float;

  -- DESCRIPTION :
  --   Returns the radius of the complex number.

  function Angle ( c : Complex_Number ) return double_float;

  -- DESCRIPTION :
  --   Returns the angle of the complex number.

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the ith root of the equation x^n - c = 0.

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number;
  function Polar_Exponentiation
             ( x : Complex_Number; e : double_float ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x.

  function Polar_Exponentiation_of_Unit
             ( x : Complex_Number; e : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x, where |x| = 1.

end Standard_Complex_Numbers_Polar;
