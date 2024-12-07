with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;

package Test_Double_Integers is

-- DESCRIPTION :
--   Tests double integer arithmetic.

  procedure Random_Double_Integer ( hi,lo : out integer64 );

  -- DESCRIPTION :
  --   Generates a random double integer, returning in hi and lo
  --   the positive high and low word.

  function Value ( hi,lo : integer64; verbose : boolean := true )
                 return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the double integer given by the
  --   high and low words in hi and lo.
  --   If verbose, then the words of the 64-bit doubles are shown.

  function Value ( hihi,lohi,hilo,lolo : integer64;
                   verbose : boolean := true ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the quadruple integer given by the
  --   words in hihi, lohi, hilo, and lolo.
  --   If verbose, then the words of the 64-bit doubles are shown.

  procedure Test_Double_Sum;

  -- DESCRIPTION :
  --   Generates two positive double integer numbers and makes their sum
  --   using double integer arithmetic.

  procedure Test_Product;

  -- DESCRIPTION :
  --   Generates two positive integer numbers and makes their product
  --   using double integer arithmetic.

  procedure Test_Double_Product;

  -- DESCRIPTION :
  --   Generates two positive double integer numbers and makes their product
  --   using double integer arithmetic.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Double_Integers;
