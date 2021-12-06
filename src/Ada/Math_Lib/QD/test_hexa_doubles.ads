with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Test_Hexa_Doubles is

-- DESCRIPTION :
--   Test operations on hexa double numbers.

  function random return hexa_double;

  -- DESCRIPTION :
  --   Returns a random octo double number from adding
  --   random double numbers in [-1,+1].

  procedure Write ( x : in hexa_double );

  -- DESCRIPTION :
  --   Writes all parts of x, one part per line.

  procedure Test_Add_and_Subtract;

  -- DESCRIPTION :
  --   Generates two random numbers, adds and subtracts.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Generates two random numbers, multiplies and divides.

  procedure Test_Read;

  -- DESCRIPTION :
  --   Reads a 256-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  procedure Log10log2exp1_doubles;

  -- DESCRIPTION :
  --   Shows the doubles for a 256-digit approximation of log(10),
  --   log(2), and exp(1).

  procedure inverse_factorials;

  -- DESCRIPTION :
  --   Prints the doubles for the inverse factorials needed in
  --   the hexa double approximation for the exp() function.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a hexa double.

  procedure Test_sqrt2;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in hexa double arithmetic.

  procedure Test_hd_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in hexa double precision,
  --   when printed with precision equal to 255.

  procedure Log_exp_of_Pi;

  -- DESCRIPTION :
  --   Tests whether log(exp(pi)) = pi = exp(log(pi)).

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Hexa_Doubles;
