with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Vectored_Octo_Doubles is

-- DESCRIPTION :
--   Test on the vectorized octo double arithmetic.

  procedure Test_Real_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random octo double numbers
  --   of dimension dim and tests their inner product.

  procedure Test_Complex_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random complex octo double numbers
  --   of dimension dim and tests their inner product.

  procedure Test_Balanced_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Tests the product of balanced quarter octo double vectors.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and runs the tests.

end Test_Vectored_Octo_Doubles;
