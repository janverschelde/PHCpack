with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;

package Multprec_Complex_Numbers_Polar is

-- DESCRIPTION :
--   Offers a polar view on the multi-precision complex numbers.

  function Radius ( c : Complex_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the radius of the complex number.

  function Angle ( c : Complex_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the angle of the complex number.

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the ith root of the equation x^n - c = 0.
  --   The size of the fraction of the real part of c determines the
  --   accuracy (i.e.: the size) of the number on return.

end Multprec_Complex_Numbers_Polar;
