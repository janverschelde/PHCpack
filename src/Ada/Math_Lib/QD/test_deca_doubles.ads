with Deca_Double_Numbers;                use Deca_Double_Numbers;

package Test_Deca_Doubles is

-- DESCRIPTION :
--   Test operations on deca double numbers.

  procedure write ( x : deca_double );

  -- DESCRIPTION :
  --   Writes all parts of x, one part per line.

  function random return deca_double;

  -- DESCRIPTION :
  --   Returns a random deca double number from adding
  --   random double numbers in [-1,+1].

  procedure Test_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Generates two random numbers, adds and subtracts.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Generates two random numbers, multiplies and divides.

  procedure Test_Read;

  -- DESCRIPTION :
  --   Reads a 160-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  procedure Test_Write;

  -- DESCRIPTION :
  --   Tests writing of a double.

  procedure Log10log2exp1_doubles;

  -- DESCRIPTION :
  --   Shows the doubles for a 160-digit approximation 
  --   of log(10), log(2), and exp(1).

  procedure inverse_factorials;

  -- DESCRIPTION :
  --   Prints the doubles for the inverse factorials needed in
  --   the deca double approximation for the exp() function.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a deca double.

  procedure Test_sqrt2;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in deca double arithmetic.

  procedure Test_da_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in deca double precision,
  --   when printed with precision equal to 159.

  procedure Write_Pi;

  -- DESCRIPTION :
  --   Writes the deca double expansion for Pi and multiples,
  --   as needed for the deca double constants.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Deca_Doubles;
