with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with PentDobl_Complex_Vector_Series;

package Test_PentDobl_Linearization is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series,
--   in penta double precision.

  procedure Write_Difference
              ( x,y : in PentDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Writes the max norm of the difference of each coefficient vector
  --   between x and y.  At the end, writes the largest max norm, as an
  --   upper bound on the error.

  -- REQUIRED : x.deg = y.deg >= 0.

  procedure PentDobl_Test ( n,m,d : in integer32 );

  -- DESCRIPTION :
  --   Generates an n-by-m matrix of series of degree d,
  --   with complex coefficients in double double precision.
  --   Converts an n-by-m matrix of series of degree d with penta
  --   double precision complex coefficients into a matrix series.

  procedure PentDobl_Timing ( n,m,d,f : in integer32 );

  -- DESCRIPTION :
  --   Generates a random problem and solves it f times
  --   by LU in case n = m or with QR if n > m,
  --   in penta double precision.

  -- ON ENTRY :
  --   n        number of equations, number of rows of the matrices;
  --   m        number of variables, number of columns of the matrices;
  --   d        degree of the series;
  --   f        frequency of tests.

  procedure PentDobl_Coefficient_Test ( n,m,d : in integer32 );

  -- DESCRIPTION :
  --   Generates an n-by-m matrix of series of degree d,
  --   with complex coefficients in penta double precision.
  --   Converts an n-by-m matrix of series of degree d with penta
  --   double precision complex coefficients into a matrix series.
  --   The linearization tested operates directly on coefficient vectors
  --   of power series, without data type encapsulation.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension of the linear system,
  --   the degrees of the series in the system,
  --   and the frequency of tests for timing.

end Test_PentDobl_Linearization;
