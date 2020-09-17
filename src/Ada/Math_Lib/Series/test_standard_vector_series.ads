with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Series_Vectors;

package Test_Standard_Vector_Series is

-- DESCRIPTION :
--   Tests vectors of truncated power series in double precisino.

  procedure Write ( v : in Standard_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the components of the vector to standard output.
 
  procedure Standard_Test_Norm 
              ( v : in Standard_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests the normalization of the vector v.

  procedure Standard_Test ( dim,degree : in integer32 );

  -- DESCRIPTION :
  --   Generates a random vector of range 1..dim,
  --   with series of the given degree, in standard double precision.
  --   Tests the computation of the norm and the normalization.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a dimension, and degree, and runs a test.

end Test_Standard_Vector_Series;
