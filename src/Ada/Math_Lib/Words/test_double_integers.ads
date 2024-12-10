with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;

package Test_Double_Integers is

-- DESCRIPTION :
--   Tests double integer arithmetic.

  procedure Random_Double52_Integer ( hi,lo : out integer64 );

  -- DESCRIPTION :
  --   Generates a random double integer, returning in hi and lo
  --   the positive high and low word, both less than 2^52.

  procedure Random_Double64_Integer ( hi,lo : out integer64 );

  -- DESCRIPTION :
  --   Generates a random double integer, returning in hi and lo
  --   the positive high and low word, both less than 2^62.

  function Value52 ( hi,lo : integer64 ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the double integer given by the
  --   high and low words in hi and lo, in base 2^52.

  function Value60 ( hi,lo : integer64; verbose : boolean := true )
                   return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the double integer given by the
  --   high and low words in hi and lo, in base 2^60.
  --   If verbose, then the words of the 64-bit doubles are shown.

  function Value60 ( hihi,lohi,hilo,lolo : integer64;
                     verbose : boolean := true ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the quadruple integer given by the
  --   words in hihi, lohi, hilo, and lolo, in base 2^60.
  --   If verbose, then the words of the 64-bit doubles are shown.

  procedure Test_Double52_Sum;

  -- DESCRIPTION :
  --   Generates two positive double integer numbers and makes their sum
  --   using double integer arithmetic, with base 2^52.

  procedure Test_Double60_Sum;

  -- DESCRIPTION :
  --   Generates two positive double integer numbers and makes their sum
  --   using double integer arithmetic, with base 2^60.

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
