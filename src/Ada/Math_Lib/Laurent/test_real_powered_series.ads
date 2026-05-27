with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Real_Powered_Series is

-- DESCRIPTION :
--   Tests the input/output procedures for real powered series
--   and the real power series arithmetic.

  procedure Write ( cff : in Standard_Complex_Vectors.Vector;
                    pwt : in Standard_Floating_Vectors.Vector;
                    t : in character := 't' );

  -- DESCRIPTION :
  --   Writes the series with coefficients in cff and the powers in pwt,
  --   using the character t as the symbol for the parameter.

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

  procedure Test_Series_IO ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random power series of the given size
  --   to test the input/output procedures.

  procedure Test_Random_Series ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random real power series of the given size.

  procedure Test_Addition ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their addition and the subtraction.

  procedure Test_Multiplication ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their multiplication.

  procedure Test_Inverse ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random series of the given size,
  --   computes the inverse and multiplies the random series
  --   with the inverse to verify the correctness.

  procedure Test_Division ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their division.

  procedure Test_Series_IO;

  -- DESCRITION :
  --   Main test on the series input and output procedures.

  procedure Test_Series_Arithmetic;

  -- DESCRIPTION :
  --   Tests the arithmetic with real power series.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Real_Powered_Series;
