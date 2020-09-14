package Test_PentDobl_Complex is

-- DESCRIPTION :
--   Tests the operations on complex numbers in penta double precision.

  procedure Test_io;

  -- DESCRIPTION :
  --   Prompts for a complex number and writes the number.

  procedure Test_Addition_and_Subtraction;

  -- DESCRIPTION :
  --   Tests x + y - x for randomly generated complex penta doubles.

  procedure Test_Multiplication_and_Division;

  -- DESCRIPTION :
  --   Tests x * y / x for randomly generated complex penta doubles.

  procedure Test_Random;

  -- DESCRIPTION :
  --   Generates a random complex number on the unit circle
  --   and computes its radius to verify it is on the unit circle.

  procedure Test_Roots;

  -- DESCRIPTION :
  --   Computes the roots of x^d - c = 0,
  --   prompting for the degree d and complex number c.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_PentDobl_Complex;
