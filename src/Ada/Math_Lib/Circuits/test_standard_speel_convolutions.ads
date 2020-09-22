with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;

package Test_Standard_Speel_Convolutions is

-- DESCRIPTION :
--   Tests a vectorized version of the computation of Speelpenning products
--   for truncated power series in standard double precision.

  procedure Standard_Test ( dim,deg : in integer32 );

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg,
  --   and runs a test in standard double precision.

  procedure Standard_Indexed_Test
              ( dim,deg,nz : in integer32;
                xp : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg, runs a test
  --   for a random exponent vector xp of zeros and ones,
  --   where the number of nonzeros nz > 2.

  procedure Indexed_Test ( dim,deg : in integer32 );

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg, and tests 
  --   for a random exponent vector of zeros and ones.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the degree, the dimension,
  --   and the setup for a test.

end Test_Standard_Speel_Convolutions;
