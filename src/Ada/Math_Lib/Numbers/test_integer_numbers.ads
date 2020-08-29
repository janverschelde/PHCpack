with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer_Numbers;

package Test_Integer_Numbers is

-- DESCRIPTION :
--   Tests operations with multiprecision integer numbers.

  procedure Test_Creation;

  -- DESCRIPTION :
  --   Tests the making of an integer number.

  procedure Test64_Creation;

  -- DESCRIPTION :
  --   Tests the making of an integer number with 64-bit words.

  procedure Test_io;

  -- DESCRIPTION :
  --   Tests the input and output.

  procedure Test64_io;

  -- DESCRIPTION :
  --   Tests the input and output of integer numbers with 64-bit words.

  procedure Test_Sign ( i : in Multprec_Integer_Numbers.Integer_Number );

  -- DESCRIPTION :
  --   Applies the operations to determine the sign of the number i.

  procedure Test_Compare
              ( i1,i2 : in Multprec_Integer_Numbers.Integer_Number );

  -- DESCRIPTION :
  --   Compares the numbers i1 and i2.
  --   Writes the result of the comparison to screen.

  procedure Zero_Test ( i : Multprec_Integer_Numbers.Integer_Number );

  -- DESCRIPTION :
  --   Tests whether the number i is zero or not.

  procedure Test_Comparison;

  -- DESCRIPTION :
  --   Tests comparison and copying operations.

  procedure Test_Addition;

  -- DESCRIPTION :
  --   Tests the addition of user given integer numbers.

  procedure Test64_Addition;

  -- DESCRIPTION :
  --   Prompts for integers and test the addition, with 64-bit words.

  function Mult_by_Add
             ( i1,i2 : Multprec_Integer_Numbers.Integer_Number )
             return Multprec_Integer_Numbers.Integer_Number;

  -- DESCRIPTION :
  --   Does the multiplication by adding up n1 to itself as many times
  --   as the number i2.  Only to be used as test of course.
  --   This can be quite time consuming as i2 gets large.

  procedure Test_Multiplication;

  -- DESCRIPTION :
  --   Tests the multiplication of two integer numbers.

  procedure Test_Exponentiation;

  -- DESCRIPTION :
  --   Tests the exponentiation.

  procedure Test64_Exponentiation;

  -- DESCRIPTION :
  --   Tests the exponentiation of integer numbers with 64-bit words.

  procedure Test_Subtraction;

  -- DESCRIPTION :
  --   Tests the subtraction of integer numbers.

  procedure Test64_Subtraction;

  -- DESCRIPTION :
  --   Tests the subtraction of integer number with 64-bit words.

  procedure Test_Division;

  -- DESCRIPTION :
  --   Tests the division of integer numbers.

  procedure Test64_Division;

  -- DESCRIPTION :
  --   Tests the division of integer numbers with 64-bit words.

  procedure Random_Addition_and_Subtraction ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?
  --   The size of the numbers is given in sz1 and sz2.

  procedure Random64_Addition_and_Subtraction ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?
  --   The size of the number is given in sz1 and sz2.

  procedure Additions_and_Subtractions_on_Randoms ( sixtyfour : in boolean );

  -- DESCRIPTION :
  --   Generates a number of random integers and performs repeated
  --   additions and subtractions with checks on consistencies.
  --   If sixtyfour, then 64-bit words are used.

  procedure Random_Multiplication_and_Division ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Four tests are performed :
  --   1) n1*n2/n2 = n1, with "*" and "/".
  --   2) Mul(n1,n2) is the same as n1 := n1*n2 ?
  --   3) Div(n1*n2,n1) leads to n2 ?
  --   4) n1 = (n1/n2)*n2 + Rmd(n1,n2) ?
  --   5) Div(n1,n2,q,r) satisfies n1 = q*n2 + r ?
  --   The size of the numbers is given in sz1 and sz2.

  procedure Random64_Multiplication_and_Division ( sz1,sz2 : in natural32 );

  -- DESCRIPTION :
  --   Four tests are performed :
  --   1) n1*n2/n2 = n1, with "*" and "/".
  --   2) Mul(n1,n2) is the same as n1 := n1*n2 ?
  --   3) Div(n1*n2,n1) leads to n2 ?
  --   4) n1 = (n1/n2)*n2 + Rmd(n1,n2) ?
  --   5) Div(n1,n2,q,r) satisfies n1 = q*n2 + r ?
  --   The size of the numbers is given in sz1 and sz2.

  procedure Multiplications_and_Divisions_on_Randoms ( six4 : in boolean );

  -- DESCRIPTION :
  --   Generates a number of random integers and performs repeated
  --   multiplications and divisions with checks on consistencies.
  --   If six4, then 64-bit words are used.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Integer_Numbers;
