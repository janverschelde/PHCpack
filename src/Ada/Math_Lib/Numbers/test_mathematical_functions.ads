with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package Test_Mathematical_Functions is

-- DESCRIPTION :
--   Tests some mathematical functions.

  function C_COS ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Calls the C cosine function.

  function C_SIN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Calls the C sine function.

  procedure Read ( f : in out Floating_Number; name : in string );

  -- DESCRIPTION :
  --   Prompts the user to read a floating number,
  --   using the name in the prompt.

  procedure Test_Standard_COS_and_SIN;

  -- DESCRIPTION :
  --   Simple test on standard sine and cosine,
  --   checks if cos^2 + sin^2 holds for a numerical input.

  procedure Test_Standard_SQRT;

  -- DESCRIPTION
  --   Test on y**2 - x = 0, where y is SQRT(x).

  procedure Test_Multprec_SQRT
              ( x : in Floating_Number; output : in boolean );

  -- DESCRIPTION
  --   Test on y**2 - x = 0, where y is SQRT(x).
  --   If output, then y and y**2 will be shown.

  procedure Interactive_Test_Multprec_SQRT;

  -- DESCRIPTION :
  --   Prompts for various numbers to test the SQRT on.

  procedure Random_Test_Multprec_SQRT;

  -- DESCRIPTION :
  --   Generates random numbers to compute SQRT of.

  procedure ln_one_plus_x ( x,eps : in Floating_Number;
                            nit : out natural32; max : in natural32;
                            y : out Floating_Number; fail : out boolean );

  -- DESCRIPTION :
  --   Test on the series for the natural logarithm of 1 + x.

  procedure Test_Series;

  -- DESCRIPTION :
  --   Checks teh convergence of the series ln(1+x).

  function Read_Base return character;

  -- DESCRIPTION :
  --   Interactive reading of the base of the logarithm.
  --   Prompts the user to enter e, 2, or 10 and returns one
  --   of the following characters: 'e', '2', or 'A'.

  procedure Test_Multprec_Logarithm;

  -- DESCRIPTION :
  --   Tests the multiprecision natural logarithm
  --   and the 2-logarithm and the 10-logarithm.

  procedure Test_Multprec_Exponential;

  -- DESCRIPTION :
  --   Tests the multiprecision exponential
  --   and the exponentiation with bases 2 and 10.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Mathematical_Functions;
