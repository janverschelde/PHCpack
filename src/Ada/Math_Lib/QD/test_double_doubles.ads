with Double_Double_Numbers;              use Double_Double_Numbers;

package Test_Double_Doubles is

-- DESCRIPTION :
--   Test operations on double doubles.

  procedure character_arithmetic;

  -- DESCRIPTION :
  --   Experiment on computing with characters in Ada,
  --   as needed for the translation of the C code.

  procedure Write ( d : in double_double );
 
  -- DESCRIPTION :
  --   Very basic output of a double double d.

  procedure Basic_Test;

  -- DESCRIPTION :
  --   Makes a double double from an integer given by the user
  --   and tests the from C imported ldexp function.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a double double.

  procedure Add_Sub_of_Pi_e;

  -- DESCRIPTION :
  --   Test on Pi + e - Pi.

  procedure Div_sqr_of_Pi;

  -- DESCRIPTION :
  --   Divides pi^2 of pi.

  procedure Log_exp_of_Pi;

  -- DESCRIPTION :
  --   Tests whether log(exp(pi)) = pi.

  procedure my_sqrt;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in double double arithmetic.

  procedure Test_Random;

  -- DESCRIPTIN :
  --   Generates and shows a random number.

  procedure Test_dd_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in double double precision,
  --   when printed with precision equal to 31.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Double_Doubles;
