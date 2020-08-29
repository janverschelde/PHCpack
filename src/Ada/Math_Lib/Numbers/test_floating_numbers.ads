with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Floating64_Numbers;

package Test_Floating_Numbers is

-- DESCRIPTION :
--   Test operations on multiprecision floating-point numbers.

  procedure Read ( f : in out Multprec_Floating_Numbers.Floating_Number;
                   name : in string );

  -- DESCRIPTION :
  --   Prompts for a floating-point number f,
  --   using the name in the prompt.

  procedure Read ( f : in out Multprec_Floating64_Numbers.Floating_Number;
                   name : in string );

 -- DESCRIPTION :
 --   Prompts for a floating-point number f with 64-bit words,
 --   using the name in the prompt.

  procedure Formatted_Output
              ( f : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Reads the format parameters and writes the floating-point number f
  --   according to those parameters.

  procedure Formatted64_Output
              ( f : in Multprec_Floating64_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Reads the format parameters and writes the floating-point number f
  --   according to those parameters.

  procedure Test_io;

  -- DESCRIPTION :
  --   Reads and writes a floating-point number.

  procedure Test64_io;

  -- DESCRIPTION :
  --   Reads and writes a floating-point number, with 64-bit words.

  procedure Test_Creation;

  -- DESCRIPTION :
  --   Tests the making of a multiprecision floating-point number.

  procedure Test64_Creation;

  -- DESCRIPTION :
  --   Tests the making of a multiprecision floating-point number
  --   with 64-bit words.

  procedure Test_Compare
              ( f1,f2 : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Compares f1 and f2 and writes the outcome to screen.

  procedure Test64_Compare
              ( f1,f2 : in Multprec_Floating64_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Compares f1 and f2 and writes the outcome to screen.

  procedure Zero_Test ( f : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Tests if f is zero.

  procedure Zero_Test ( f : in Multprec_Floating64_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Tests if f is zero.

  procedure Test_Comparison;

  -- DESCRIPTION :
  --   Prompts for two multiprecision floating-point numbers
  --   and tests the comparison operations.

  procedure Test64_Comparison;

  -- DESCRIPTION :
  --   Prompts for two multiprecision floating-point numbers with 64-bit words
  --   and tests the comparison operations.

  procedure Test_Size;

  -- DESCRIPTION :
  --   Tests the operations to truncate, round, and expand the precision.

  procedure Test_Addition;

  -- DESCRIPTION :
  --   Tests the addition of multiprecision floating-point numbers.

  procedure Test64_Addition;

  -- DESCRIPTION :
  --   Prompts for multiprecision floating-point numbers with 64-bit words
  --   and tests the addition.

  procedure Test_Subtraction;

  -- DESCRIPTION :
  --   Tests the subtraction of multiprecision floating-point numbers.

  procedure Test64_Subtraction;

  -- DESCRIPTION :
  --   Prompts for multiprecision floating-point numbers with 64-bit words
  --   and tests the subtraction.

  procedure Test_Multiplication;

  -- DESCRIPTION :
  --   Tests the multiplication of multiprecision floating-point numbers.

  procedure Test_Exponentiation;

  -- DESCRIPTION :
  --   Tests the exponentiation of multiprecision floating-point numbers.
  
  procedure Test_Division;

  -- DESCRIPTION :
  --   Tests the division of multiprecision floating-point numbers.

  function Min ( a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  function Max ( a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  function Acceptable_Tolerance
             ( sz1,sz2 : natural32 )
             return Multprec_Floating_Numbers.Floating_Number;

  -- DESCRIPTION :
  --   Returns an acceptable tolerance for the error of (n1+n2)-n1.
  --   If sz = min(sz1,sz2) + 1, then tol = Radix**(-sz*Exponent+2).

  procedure Report_Bug ( n1,n2 : Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Prints the numbers n1 and n2 for use in a bug report.

  procedure Compare_with_Tolerance
              ( n1,n2,resi : in Multprec_Floating_Numbers.Floating_Number; 
                acctol : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   If the residual resi on the result of an operation with n1 and n2
  --   is larger than the acceptable tolerance in acctol,
  --   then n1 and n2 are reported as bug cases.

  procedure Random_Addition_and_Subtraction
              ( sz1,sz2 : in natural32; low,upp : in integer32;
                acctol : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?
  --   The size of the numbers is in sz1 and sz2.
  --   The number low and up are bounds for the random numbers.
  --   The acctol is the acceptable tolerance on the result.

  procedure Additions_and_Subtractions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random floats and performs repeated
  --   additions and subtractions with checks on consistencies.

  procedure Random_Multiplication_and_Division
               ( sz1,sz2 : in natural32; low,upp : in integer32;
                 acctol : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   Three tests are performed :
  --   1) n1*n2/n2 = n1, with "*" and "/".
  --   2) Mul(n1,n2) is the same as n1 := n1*n2 ?
  --   3) Div(n1*n2,n1) leads to n2 ?
  --   The size of the numbers is in sz1 and sz2.
  --   The number low and up are bounds for the random numbers.
  --   The acctol is the acceptable tolerance on the result.

  procedure Multiplications_and_Divisions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random floats and performs repeated
  --   multiplications and divisions with checks on consistencies.

  procedure Multiplication_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and multiplying them.
  --   Timing is reported at the end of the cycle.

  procedure Multiplication64_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and multiplying them.
  --   Computations are done on 64-bit words.
  --   Timing is reported at the end of the cycle.

  procedure Division_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and dividing them.
  --   Timing is reported at the end of the cycle.

  procedure Division64_Timer;

  -- DESCRIPTION :
  --   Starts a cycle of generating random numbers and dividing them.
  --   Computations are done on 64-bit words.
  --   Timing is reported at the end of the cycle.

  procedure Nearest_Integer;

  -- DESCRIPTION :
  --   Interactive test on the computation of the nearest integer.

  procedure Nearest_Integer64;

  -- DESCRIPTION :
  --   Interactive test on the computation of the nearest integer,
  --   for multiprecision floating-point numbers with 64-bit words.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Floating_Numbers;
