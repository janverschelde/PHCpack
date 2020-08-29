with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;

package Test_DoblDobl_Complex is

-- DESCRIPTION :
--   Tests complex arithmetic in double double precision.

  procedure Basic_Test;

  -- DESCRIPTION :
  --   Tests input and output of a complex number.

  procedure Test_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Tests x + y - x for randomly generated complex double doubles.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Tests x * y / x for randomly generated complex double doubles.

  procedure Prompt_Complex_Number ( c : out Complex_Number );

  -- DESCRIPTION :
  --   Prompts for the real and imaginary part of a complex number
  --   returned in c.

  procedure Test_Roots;

  -- DESCRIPTION :
  --   Computes the roots of x^d - c = 0,
  --   prompting for the degree d and complex number c.

  procedure Test_Random;

  -- DESCRIPTION :
  --   Generates a random number and shows it.
 
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_DoblDobl_Complex;
