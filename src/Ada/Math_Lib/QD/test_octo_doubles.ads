with Octo_Double_Numbers;                use Octo_Double_Numbers;

package Test_Octo_Doubles is

-- DESCRIPTION :
--   Test operations on octo double numbers.

  function random return octo_double;

  -- DESCRIPTION :
  --   Returns a random octo double number from adding
  --   random double numbers in [-1,+1].

  procedure Write ( x : in octo_double );

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
  --   Reads a 128-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  procedure Log10log2exp1_doubles;

  -- DESCRIPTION :
  --   Shows the doubles for a 128-digit approximation of log(10),
  --   log(2), and exp(1).

  procedure inverse_factorials;

  -- DESCRIPTION :
  --   Prints the doubles for the inverse factorials needed in
  --   the octo double approximation for the exp() function.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a octo double.

  procedure Test_sqrt2;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in octo double arithmetic.

  procedure Test_od_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in octo double precision,
  --   when printed with precision equal to 127.

  procedure Log_exp_of_Pi;

  -- DESCRIPTION :
  --   Tests whether log(exp(pi)) = pi = exp(log(pi)).

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Octo_Doubles;
