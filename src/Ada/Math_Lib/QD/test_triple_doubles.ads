with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

package Test_Triple_Doubles is

-- DESCRIPTION :
--   Tests operations on triple double numbers.

  function to_triple_double ( x : quad_double ) return triple_double;

  -- DESCRIPTION :
  --   Copies the first three words of x into the result.

  function create ( x : triple_double ) return quad_double;

  -- DESCRIPTION :
  --   Copies the words of x into the first three words of the result.

  procedure Test_Basic_Arithmetic;

  -- DESCRIPTION :
  --   Generates random numbers and compares with quad double arithmetic.
  --   The sin and cos functions are applied to have four nonzero words.

  procedure Test_Read;

  -- DESCRIPTION :
  --   Reads a 50-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(50)
  --   >>> s2
  --   1.4142135623730950488016887242096980785696718753769
  --   >>> s2*s2 - 2
  --   -2.6727647100921956461405364671514818788151968801050e-51
    
  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a triple double.

  procedure Test_sqrt2;

  -- DESCRIPTION :
  --   Computes the square root of 2 using Newton's method
  --   in triple double arithmetic.

  procedure Test_td_eps;

  -- DESCRIPTION :
  --   Tests on the smallest number which still makes a difference
  --   when added to one in triple double precision,
  --   when printed with precision equal to 47.

  procedure Log_exp_of_Pi;

  -- DESCRIPTION :
  --   Tests whether log(exp(pi)) = pi = exp(log(pi)).

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user to select a test
  --   and then runs the test.

end Test_Triple_Doubles;
