with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Vectored_Hexa_Doubles is

-- DESCRIPTION :
--   Test on the vectorized hexa double arithmetic.

  procedure Test_Balanced_Product ( dim,freq : in integer32 );

  -- DESCRIPTION :
  --   Tests the product of balanced quarter hexa double vectors,
  --   on vectors of dimension dim, computing it as many times
  --   as the value of freq.

  procedure Wall_Time_Test;

  -- DESCRIPTION :
  --   Runs a test without prompting for input,
  --   suitable to measure the wall clock time.

  procedure Wall_Time_Parallel_Test ( nt : in integer32 );

  -- DESCRIPTION :
  --   Runs the Wall_Time_Test with a number of threads equal to nt.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and runs the tests.

end Test_Vectored_Hexa_Doubles;
