package Test_Integer_Vectors is

-- DESCRIPTION :
--   Tests vectors of 32-bit, 64-bit, and multiprecision integer numbers.

  procedure Test_Standard_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector
  --   of 32-bit integer numbers.  Writes the vector.

  procedure Test_Standard64_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector
  --   of 64-bit integer numbers.  Writes the vector.

  procedure Test_Standard_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector of vectors
  --   of 32-bit integer numbers.  Writes the vector.

  procedure Test_Standard64_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector of vectors
  --   of 64-bit integer numbers.  Writes the vector.

  procedure Test_Multprec_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for a vector of multiprecision
  --   integer numbers.  Writes the vector.

  procedure Test_Multprec_VecVecs_io;

  -- DESCRIPTION :
  --   Prompts for a dimension and then for vectors of vectors
  --   of multiprecision integer numbers.  Writes the vector.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Integer_Vectors;
