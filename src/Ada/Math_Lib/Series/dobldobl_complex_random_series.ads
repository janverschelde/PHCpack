with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Series;            use DoblDobl_Complex_Series;

package DoblDobl_Complex_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given degree,
--   with coefficients in double double precision.

  function Random_Series ( degree : integer32 ) return Series;
  function Random_Series ( degree : integer32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree, with random coefficients,
  --   on the unit circle on the complex plane.

end DoblDobl_Complex_Random_Series;
