with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;

package Test_Complex_Numbers is

-- DESCRIPTION :
--   Tests operations on standard double and multiprecision complex numbers.

  procedure Test_Standard_io;

  -- DESCRIPTION :
  --   Tests the input and output of complex numbers
  --   in standard double precision.

  procedure Test_Multprec_io;

  -- DESCRIPTION :
  --   Tests the input and output of complex multiprecision numbers.

  procedure Test_Standard_Roots;

  -- DESCRIPTION :
  --   Solves x^d - c = 0 after prompting
  --   for a degree d and complex number c,
  --   in standard double precision.

  procedure Test_Multprec_Roots;

  -- DESCRIPTION :
  --   Solves x^d - c = 0 after prompting
  --   for a degree d, a complex number c,
  --   and the size of the multiprecision numbers.

  procedure Standard_Random_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?

  procedure Standard_Additions_and_Subtractions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random floats and performs repeated
  --   additions and subtractions with checks on consistencies.

  function Min ( a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  function Acceptable_Tolerance
              ( sz1,sz2 : natural32 )
              return Multprec_Floating_Numbers.Floating_Number;

  -- DESCRIPTION :
  --   Returns an acceptable tolerance for the error of (n1+n2)-n1.
  --   If sz = min(sz1,sz2) + 1, then tol = Radix**(-sz*Exponent+2).

  procedure Report_Bug ( n1,n2 : Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Reports the numbers n1 and n2 as bug cases.

  procedure Compare_with_Tolerance 
              ( n1,n2,resi : in Multprec_Complex_Numbers.Complex_Number; 
                acctol : in Multprec_Floating_Numbers.Floating_Number );

  -- DESCRIPTION :
  --   If the residual resi on the result of an operation with n1 and n2
  --   is larger than the acceptable tolerance in acctol,
  --   then n1 and n2 are reported as bug cases.

  procedure Multprec_Random_Addition_and_Subtraction
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

  procedure Multprec_Additions_and_Subtractions_on_Randoms;

  -- DESCRIPTION :
  --   Generates a number of random floats and performs repeated
  --   additions and subtractions with checks on consistencies.

  procedure Simul_Div ( x : in out Standard_Complex_Numbers.Complex_Number;
                        y : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Simulated division of standard double numbers,
  --   with print statements to track the order of the operations.

  procedure Simul_Div ( x : in out Multprec_Complex_Numbers.Complex_Number;
                        y : in Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Simulated division of multiprecision numbers,
  --   with print statements to track the order of the operations.

  procedure Interactive_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Prompts for numbers and tests the multiplication and division.

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

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Complex_Numbers;
