with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Series_Vectors;

package Test_DoblDobl_Vector_Series is

-- DESCRIPTION :
--   Tests vectors of truncated power series in double double precision.

  procedure Write ( v : in DoblDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the components of the vector to standard output.
 
  procedure DoblDobl_Test_Norm 
              ( v : in DoblDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests the normalization of the vector v.

  procedure DoblDobl_Test ( dim,degree : in integer32 );

  -- DESCRIPTION :
  --   Generates a random vector of range 1..dim,
  --   with series of the given degree, in double double precision.
  --   Tests the computation of the norm and the normalization.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for a dimension, and degree, and runs a test.

end Test_DoblDobl_Vector_Series;
