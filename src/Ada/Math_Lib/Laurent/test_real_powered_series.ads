with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Real_Powered_Series is

-- DESCRIPTION :
--   Tests the input and output procedures for real powered series.

  procedure Random_Series
              ( size : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a random series with positive powers of the
  --   given size and sorts the powers.

  procedure Test_String_Series ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random power series of the given size
  --   to test the string conversion procedures 

  procedure Test_Random_Series ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random power series of the given size
  --   to test the input/output procedures.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Real_Powered_Series;
