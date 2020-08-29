with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;

package Test_QuadDobl_Complex is

-- DESCRIPTION :
--   Tests on complex arithmetic in double double precision.

  procedure Basic_Test;

  -- DESCRIPTION :
  --   Tests the input and output.

  procedure Test_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Tests x + y - x for randomly generated complex quad doubles.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Tests x * y / x for randomly generated complex double doubles.

  procedure Prompt_Complex_Number ( c : out Complex_Number );

  -- DESCRIPTION :
  --   Prompts for the real and imaginary part of a complex number
  --   returned in c.

  procedure Test_Roots;

  -- DESCRIPTION :
  --   Solves x^d - c = 0, after prompting for a degree d
  --   and a complex number c.

  procedure Test_Random;

  -- DESCRIPTION :
  --   Generates a random number and shows it.
 
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_QuadDobl_Complex;
