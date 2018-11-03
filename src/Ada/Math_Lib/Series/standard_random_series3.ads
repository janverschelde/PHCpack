with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Series;            use Standard_Complex_Series;

package Standard_Random_Series3 is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given degree,
--   with coefficients in standard double precision.

  function Random_Series ( degree : integer32 ) return Series;
  function Random_Series ( degree : integer32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree, with random coefficients,
  --   on the unit circle on the complex plane.

end Standard_Random_Series3;
