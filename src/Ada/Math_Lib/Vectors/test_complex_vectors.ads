package Test_Complex_Vectors is

-- DESCRIPTION :
--   Tests vectors of standard double and multiprecision complex numbers.

  procedure Test_Standard_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for the dimension and then for a complex vector
  --   of that dimension in standard double precision.

  procedure Test_Multprec_Vectors_io;

  -- DESCRIPTION :
  --   Prompts for the dimension and then for a complex vector
  --   of that dimension in multiprecision.

  procedure Test_Random_Vectors;

  -- DESCRIPTION :
  --   Prompts for the dimension and the magnitude of the numbers
  --   and then shows a randomly generated vector of complex numbers
  --   in standard double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Complex_Vectors;
