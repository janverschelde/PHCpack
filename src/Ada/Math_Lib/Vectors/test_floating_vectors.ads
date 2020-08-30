package Test_Floating_Vectors is

-- DESCRIPTION :
--   Tests on vectors of standard double and multiprecision floats.

  procedure Test_Standard_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for the entries
  --   of a vector of multiprecision floating-point numbers.

  procedure Test_Standard_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for the dimension and then prompts for vectors of vectors
  --   of standard double precision floating-point numbers.

  procedure Test_Multprec_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for the entries
  --   of a vector of multiprecision floating-point numbers.

  procedure Test_Multprec_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for the dimension and then prompts for vectors of vectors
  --   of multiprecision floating-point numbers.

  procedure Test_Random_Vectors;

  -- DESCRIPTION :
  --   Prompts for the dimension and the magnitude of the numbers and
  --   then shows a randomly generated vector in standard double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Floating_Vectors;
