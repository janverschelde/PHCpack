with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;

package Test_Natural_Numbers is

-- DESCRIPTION :
--   Tests operations on multiprecision natural numbers.

  procedure Test_Creation;

  -- DESCRIPTION :
  --   Tests the making of a natural number.

  procedure Test64_Creation;

  -- DESCRIPTION :
  --   Tests the making of a natural number with 64-bit words.

  procedure Test_io;

  -- DESCRIPTION :
  --   Tests input and output.

  procedure Test64_io;

  -- DESCRIPTION :
  --   Tests input and output of a natural number with 64-bit words.

  procedure Test_Compare
              ( n1,n2 : in Multprec_Natural_Numbers.Natural_Number );

  -- DESCRIPTION :
  --   Compares two natural numbers n1 and n2,
  --   writes the result of the comparison.

  procedure Test_Comparison;

  -- DESCRIPION :
  --   Prompts for two natural numbers and compares them.

  procedure Test_Addition;

  -- DESCRIPTION :
  --   Tests the addition of two natural numbers.

  procedure Test64_Addition;

  -- DESCRIPTION :
  --   Tests the addition of two natural numbers with 64-bit words.

  function Mult_by_Add ( n1,n2 : Multprec_Natural_Numbers.Natural_Number )
                       return Multprec_Natural_Numbers.Natural_Number;

  -- DESCRIPTION :
  --   Does the multiplication by adding up n1 to itself as many times
  --   as the number n2.  Only to be used as test of course.
  --   This can be quite time consuming as n2 gets large.

  procedure Test_Multiplication;

  -- DESCRIPTION :
  --   Tests the multiplication.

  procedure Test_Exponentiation;

  -- DESCRIPTION :
  --   Tests the exponentation.

  procedure Test_Subtraction;

  -- DESCRIPTION :
  --   Tests subtraction.

  procedure Divide10 ( n : in Multprec_Natural_Numbers.Natural_Number );

  -- DESCRIPTION :
  --   Checks whether the number n is divisible by 1..10.

  procedure Test_Short_Division;

  -- DESCRIPTION :
  --   This routine tests the division by a small divisor.

  procedure Test_Long_Division;

  -- DESCRIPTION :
  --   This routine tests the division by a long divisor.

  procedure Test_Division ( short : in boolean );

  -- DESCRIPTION :
  --   If short, then the divisor is a standard natural number.

  procedure Test64_Division;

  -- DESCRIPTION :
  --   Tests division with natural numbers of 64-bit words.

  procedure Random_Addition_and_Subtraction ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?

  procedure Additions_and_Subtractions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random naturals and performs repeated
  --   additions and subtractions with checks on consistencies.

  procedure Random_Multiplication_and_Division ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Four tests are performed :
  --   1) n1*n2/n2 = n1, with "*" and "/".
  --   2) Mul(n1,n2) is the same as n1 := n1*n2 ?
  --   3) Div(n1*n2,n1) leads to n2 ?
  --   4) n1 = (n1/n2)*n2 + Rmd(n1,n2) ?
  --   5) Div(n1,n2,q,r) satisfies n1 = q*n2 + r ?

  procedure Random64_Multiplication_and_Division ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Four tests are performed :
  --   1) n1*n2/n2 = n1, with "*" and "/".
  --   2) Mul(n1,n2) is the same as n1 := n1*n2 ?
  --   3) Div(n1*n2,n1) leads to n2 ?
  --   4) n1 = (n1/n2)*n2 + Rmd(n1,n2) ?
  --   5) Div(n1,n2,q,r) satisfies n1 = q*n2 + r ?

  procedure Multiplications_and_Divisions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random naturals and performs repeated
  --   multiplications and divisions with checks on consistencies.

  procedure Multiplications64_and_Divisions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random naturals and performs repeated
  --   multiplications and divisions with checks on consistencies.

  procedure Multiplication_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and multiplying them.
  --   Timing is reported at the end of the cycle.

  procedure Multiplication64_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and multiplying them.
  --   Timing is reported at the end of the cycle.

  procedure Division_Remainder_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and dividing them.
  --   Timing is reported at the end of the cycle.

  procedure Division64_Remainder_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and dividing them.
  --   Timing is reported at the end of the cycle.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Natural_Numbers;
