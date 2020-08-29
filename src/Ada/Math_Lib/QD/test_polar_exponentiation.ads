with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Test_Polar_Exponentiation is

-- DESCRIPTION :
--   Test on exponentiation of complex numbers via polar representation.

  procedure Test_Standard_DivModTwoPi;

  -- DESCRIPTION :
  --   Interactive test to compute division and remainder
  --   of any integer number by 2 times Pi, in standard double precision.

  procedure Test_DoblDobl_DivModTwoPi;

  -- DESCRIPTION :
  --   Interactive test to compute division and remainder
  --   of any integer number by 2 times Pi, in double double precision.

  procedure Test_QuadDobl_DivModTwoPi;

  -- DESCRIPTION :
  --   Interactive test to compute division and remainder
  --   of any integer number by 2 times Pi, in quad double precision.

  procedure Test_Multprec_DivModTwoPi;

  -- DESCRIPTION :
  --   Interactive test to compute division and remainder
  --   of any floating number by 2 times Pi, in arbitrary multiprecision.

  function Binary_Exponentiation 
             ( x : Complex_Number; e : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Uses binary exponentiation to compute x^e,
  --   in standard double precision.

  procedure Test_Standard_Exponentiation;

  -- DESCRIPTION :
  --   Interactive test on the exponentiation of a complex number.

  procedure Test_Multprec_Exponentiation;

  -- DESCRIPTION :
  --   Interactive test on exponentiation of a complex number
  --   with a multiprecision integer exponent.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Polar_Exponentiation;
