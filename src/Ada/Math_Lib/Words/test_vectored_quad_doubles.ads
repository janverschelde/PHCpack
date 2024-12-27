with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Vectored_Quad_Doubles is

-- DESCRIPTION :
--   Test on the vectorized quad double arithmetic.

  procedure Test_Real_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random quad double numbers
  --   of dimension dim and tests their inner product.

  procedure Test_Complex_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random complex quad double numbers
  --   of dimension dim and tests their inner product.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and runs the tests.

end Test_Vectored_Quad_Doubles;