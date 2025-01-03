with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Vectored_Double_Doubles is

-- DESCRIPTION :
--   Test on the vectorized double double arithmetic.

  procedure Test_Real_Sum ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random vector of double double numbers
  --   of dimension dim and tests their sum.

  procedure Test_Complex_Sum ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random vector of double double numbers
  --   of dimension dim and tests their sum.

  procedure Test_Real_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random double double numbers
  --   of dimension dim and tests their inner product.

  procedure Test_Complex_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two vectors of random complex double double numbers
  --   of dimension dim and tests their inner product.

  procedure Test_Real_Norm ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates one vector of random double double numbers
  --   of dimension dim and tests its squared complex norm.

  procedure Test_Complex_Norm ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates one vector of random complex double double numbers
  --   of dimension dim and tests its squared complex norm.

  procedure Test_Balanced_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Tests the product of balanced quarter double double vectors.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and runs the tests.

end Test_Vectored_Double_Doubles;
