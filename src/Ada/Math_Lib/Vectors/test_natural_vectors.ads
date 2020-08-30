package Test_Natural_Vectors is

-- DESCRIPTION :
--   Tests vectors of 32-bit natural and multiprecision natural numbers.

  procedure Test_Standard_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector
  --   of 32-bit natural numbers.

  procedure Test_Standard_Addition;

  -- DESCRIPTION :
  --   Prompts for a dimension and two vectors of 32-bit natural numbers,
  --   displays the result of the sum of the two vectors.

  procedure Test_Standard_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then prompts for vectors
  --   of vectors of 32-bit long natural numbers.

  procedure Test_Multprec_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then prompts for a vector
  --   of multiprecision natural numbers.

  procedure Test_Multprec_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then prompts for vectors
  --   of vectors of multiprecision natural numbers.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Natural_Vectors;
