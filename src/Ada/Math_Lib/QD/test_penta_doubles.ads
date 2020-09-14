with Penta_Double_Numbers;               use Penta_Double_Numbers;

package Test_Penta_Doubles is

-- DESCRIPTION :
--   Tests operations on penta doubles.

  procedure write ( x : penta_double );

  -- DESCRIPTION :
  --   Writes all parts of x, one part per line.

  function random return penta_double;

  -- DESCRIPTION :
  --   Returns a random penta double number from adding
  --   random double numbers in [-1,+1].

  procedure Test_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Generates two random numbers, adds and subtracts.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Generates two random numbers, multiplies and divides.

  procedure Test_Read;

  -- DESCRIPTION :
  --   Reads a 80-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a penta double.

  procedure Test_sqrt2;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in penta double arithmetic.

  procedure Test_pd_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in penta double precision,
  --   when printed with precision equal to 79.

  procedure Log_exp_of_Pi;

  -- DESCRIPTION :
  --   Tests whether log(exp(pi)) = pi = exp(log(pi)).

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Penta_Doubles;
