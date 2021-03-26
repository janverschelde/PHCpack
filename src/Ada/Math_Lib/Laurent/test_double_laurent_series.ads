with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Double_Laurent_Series is

-- DESCRIPTION :
--   Performs some basic tests on truncated Laurent power series,
--   defined by a leading exponent (which may be negative)
--   and a complex coefficient vector, in double precision.

  procedure Test_Multiply_Inverse_Divide ( deg : in integer32 );

  -- DESCRIPTION :
  --   Generates coefficient vectors for truncated Laurent series
  --   to degree deg, to test multiplication, inverse, and division.

  procedure Test_Add_and_Subtract ( deg : in integer32 );

  -- DESCRIPTION :
  --   Generates coefficient vectors for truncated Laurent series
  --   to degree deg, to test addition and subtraction.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the degree of truncation,
  --   and tests some basic arithmetical operations.

end Test_Double_Laurent_Series;
