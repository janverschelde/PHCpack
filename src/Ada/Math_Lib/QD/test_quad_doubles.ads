with Quad_Double_Numbers;                use Quad_Double_Numbers;

package Test_Quad_Doubles is

-- DESCRIPTION :
--   Interactive testing on the operations on quad doubles.

  procedure Write ( q : in quad_double );
 
  -- DESCRIPTION :
  --   Very basic output of a quad double q.

  procedure Basic_Test;

  -- DESCRIPTION :
  --   Makes a double double from an integer given by the user.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a quad double.

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
  --   in quad double arithmetic.

  procedure Test_Random;

  -- DESCRIPTION :
  --   Generates a random number and shows it.

  procedure Test_qd_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in quad double precision,
  --   when printed with precision equal to 63.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Quad_Doubles;
