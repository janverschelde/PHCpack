with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Greatest_Common_Divisors is

-- DESCRIPTION :
--   Interactive and random tests of gcd-computations, 
--   for 32-bit, 64-bit, and multiprecision integer numbers.

  procedure Interactive_Test_Standard_GCD;

  -- DESCRIPTION :
  --   Prompts for two 64-bit integer numbers,
  --   computes the least common multiple and 
  --   the extended greatest common divisors.

  procedure Test_Standard_GCD ( low,upp : in integer64 );

  -- DESCRIPTION :
  --   Generates two numbers between low and upp and performs the
  --   following checks :
  --    1) gcd(a,b) = d, same d as in gcd(a,b,k,l,d) ?
  --    2) Is k*a + l*b = d, with a,b,k,l,d from gcd(a,b,k,l,d) ?

  procedure Random_Test_Standard_GCD;

  -- DESCRIPTION :
  --   Prompts for the number of tests on 64-bit greatest common divisors
  --   and the lower and upper bound on the random numbers.

  procedure Interactive_Test_Multprec_GCD;

  -- DESCRIPTION :
  --   Prompts for multiprecision integer numbers,
  --   computes the least common multiple and 
  --   the extended greatest common divisors.

  procedure Test_Multprec_GCD ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Generates two numbers between low and upp and performs the
  --   following checks :
  --    1) gcd(a,b) = d, same d as in gcd(a,b,k,l,d) ?
  --    2) Is k*a + l*b = d, with a,b,k,l,d from gcd(a,b,k,l,d) ?

  procedure Random_Test_Multprec_GCD;

  -- DESCRIPTION :
  --   Prompts for the number of tests on the greatest common divisors
  --   on multiprecision integers, and on the sizes of the numbers.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Greatest_Common_Divisors;
