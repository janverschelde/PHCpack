with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Test_Random_Balanced_Quarters is

-- DESCRIPTION :
--   Test the operations on balanced quarter doubles.

  procedure Test_Thirteen_Bits;

  -- DESCRIPTION :
  --   Tests the generation of random 13-bit integers.

  procedure Write_Quarters ( x0,x1,x2,x3 : in double_float );

  -- DESCRIPTION :
  --   Writes the quarters, their sums, and their binary expansions.

  procedure Test_Random_Quarters;

  -- DESCRIPTION :
  --   Tests the generation of random balanced quarter doubles.

  procedure Test_Random_Vectors;

  -- DESCRIPTION :
  --   Tests the generation of vectors of random quarter doubles.

  procedure Test_Double_Wrapper;

  -- DESCRIPTION :
  --   Generates a random balanced quarter double and verifies
  --   whether its split in quarters is balanced.

  procedure Test_Double_Double_Wrapper;

  -- DESCRIPTION :
  --   Generates a random balanced quarter double double and verifies
  --   whether its split in quarters is balanced.

  procedure Test_Quad_Double_Wrapper;

  -- DESCRIPTION :
  --   Generates a random balanced quarter quad double and verifies
  --   whether its split in quarters is balanced.

  procedure Test_Octo_Double_Wrapper;

  -- DESCRIPTION :
  --   Generates a random balanced quarter octo double and verifies
  --   whether its split in quarters is balanced.

  procedure Test_Hexa_Double_Wrapper;

  -- DESCRIPTION :
  --   Generates a random balanced quarter hexa double and verifies
  --   whether its split in quarters is balanced.

  procedure Test_Balanced_Split;

  -- DESCRIPTION :
  --   Generates a random double and tests the splitting on a grid.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Random_Balanced_Quarters;
